/*
Description:
A simple spectrum analyzer application using a LimeSDR mini, FFTW and GNU Plot.

Author:
Fabian Török

V1.0
 */
#include "lime/LimeSuite.h"     // Required to drive LimeSDR
#include <stdlib.h>
#include <stdio.h>
#include <gnuplot_i.hpp>          // Required for Plotting
#include <time.h>
#include <math.h>
#include <fftw3.h>              // Required for DFT and FFT
#include <csvparser.h>
#include <vector>
#include <Header.h>
#include <unistd.h>
#ifdef USE_GNU_PLOT
#include "gnuPlotPipe.h"
#endif

#define MAX_TRACE 1

class Lime_SDR_mini_Rx
{
    public:
        lms_device_t* device;
        lms_stream_t streamId;
        Lime_SDR_mini_Rx()
        {
            device = NULL;
        };
        ~Lime_SDR_mini_Rx()
        {
            Close_Stream();
            LMS_Close(device);
        };


        int Setup_SDR(float_type centerFreq, float_type sampleRate, unsigned gain)
        {
            //open the first device
            lms_info_str_t list[8];
            int n;

            if ((n = LMS_GetDeviceList(list)) < 0)  //NULL can be passed to only get number of devices
                error("No devices found", device);

            if (LMS_Open(&device, list[0], NULL))
                error("Device could not be opened", device);

            if (LMS_Init(device) != 0)
                error("Device could not be initialized", device);

            //Enable RX channel
            //Channels are numbered starting at 0
            if (LMS_EnableChannel(device, LMS_CH_RX, 0, true) != 0)
                error("RX Channel could not be enabled", device);

            //Set center frequency
            if (LMS_SetLOFrequency(device, LMS_CH_RX, 0, centerFreq) != 0)
                error("Center Frequency could not be set", device);

            //Set sample rate, ask to use 2x oversampling in RF
            //This sets sampling rate for all channels
            if (LMS_SetSampleRate(device, sampleRate, 0) != 0)
                error("Sample rate could not be set", device);

            //Set Lowpass Filter Bandwidth
            LMS_SetLPFBW(device, LMS_CH_RX, 0, sampleRate*2);

            if(LMS_SetGaindB(device, LMS_CH_RX, 0, gain)!=0)
                error("Gain could not be set", device);

            return 1;
        };

        void Change_RF_Settings(float_type centerFreq = 0, float_type sampleRate = 0, unsigned gain = 100)
        {
            if(centerFreq != 0)
            {
                //Set center frequency
                if (LMS_SetLOFrequency(device, LMS_CH_RX, 0, centerFreq) != 0)
                    error("Center Frequency could not be set", device);
            }
            if(sampleRate != 0)
            {
                //Set sample rate, ask to use 2x oversampling in RF
                //This sets sampling rate for all channels
                if (LMS_SetSampleRate(device, sampleRate, 0) != 0)
                    error("Sample rate could not be set", device);
            }
            if(gain != 100)
            {
                if(LMS_SetGaindB(device, LMS_CH_RX, 0, gain)!=0)
                    error("Gain could not be set", device);
            }
        };

        int Setup_Stream()
        {
            streamId.channel = 0;               //channel number
            streamId.fifoSize = 1024 * 1024;    //fifo size in samples
            streamId.throughputVsLatency = 1.0; //optimize for max throughput
            streamId.isTx = false;              //RX channel
            streamId.dataFmt = lms_stream_t::LMS_FMT_I12;     //12-bit integers //lms_stream_t::

            //Setup Stream
            if (LMS_SetupStream(device, &streamId) != 0)
                error("Stream could not be set up", device);

            return 1;
        };
        void Start_Stream()
        {
            //Start streaming
            LMS_StartStream(&streamId);
        };
        void Stop_Stream()
        {
            //Stop streaming
            LMS_StopStream(&streamId);
        };
        void Close_Stream()
        {
            LMS_DestroyStream(device, &streamId);
        };

};

class FFT_Settings
{
    private:
        fftw_plan p;

    public:
        int DFTlength;
        fftw_complex *in;
        fftw_complex *out;

        FFT_Settings(int DFTlengthIn)
        {
            DFTlength = DFTlengthIn;
            in = (fftw_complex*) fftw_malloc(DFTlength * sizeof(fftw_complex));
            out = (fftw_complex*) fftw_malloc(DFTlength * sizeof(fftw_complex));
            p = fftw_plan_dft_1d(DFTlength, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
        };
        void Execute_FFTW()
        {
            fftw_execute(p);
        };
        void Change_FFT_Settings(int DFTlengthIn)
        {
            //Free old memory
            fftw_free(in);
            fftw_free(out);
            fftw_destroy_plan(p);

            // Set new values
            DFTlength = DFTlengthIn;
            in = (fftw_complex*) fftw_malloc(DFTlength * sizeof(fftw_complex));
            out = (fftw_complex*) fftw_malloc(DFTlength * sizeof(fftw_complex));
            p = fftw_plan_dft_1d(DFTlength, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
        };


        ~FFT_Settings()
        {
            fftw_destroy_plan(p);
            fftw_free(in);
            fftw_free(out);
        };
};


class Trace
{
    public:
    Trace(float_type centerFreqIn, float_type dispSpanIn,  float_type resBWIn, double gainIn, double refAmpIn, int FFTIn, unsigned TraceNum)
    {
        Display.centerFreq = centerFreqIn;
        Display.span = dispSpanIn;
        Display.resBW = resBWIn;
        Display.refAmp = refAmpIn;
        SDRSetting.sampleRate = resBWIn*2;
        SDRSetting.frameTime = 1/SDRSetting.sampleRate*FFTIn;
        SDRSetting.centerFreq = centerFreqIn - dispSpanIn / 2;
        SDRSetting.gain = gainIn;
        FFT.Change_FFT_Settings(FFTIn);
        YAxisBuffer = new double [int(ceil(FFTIn/2*(Display.span/Display.resBW)))];
        XAxisBuffer = new double [int(ceil(FFTIn/2*(Display.span/Display.resBW)))];
        for(int x=0; x<(int)(ceil(FFTIn/2*(Display.span/Display.resBW))); x++)
            YAxisBuffer[x] = 0;
        for(int x=0; x<int(ceil(FFTIn/2*(Display.span/Display.resBW))); x++)
            XAxisBuffer[x] = SDRSetting.centerFreq + (SDRSetting.sampleRate / FFT.DFTlength)*x;
        StreamBuffer = new int16_t[FFT.DFTlength * 2];
        SDR.Setup_SDR(SDRSetting.centerFreq, SDRSetting.sampleRate, SDRSetting.gain);
    };

    ~Trace()
    {
        free(YAxisBuffer);
        free(XAxisBuffer);
        free(StreamBuffer);
    };

    void Aqcuire_Sweep()
    {
        FILE* plotFile;
        //gnuplotHandle2 = Gnuplot("lines");
        int FFTLength = getFFTLength();
        SDR.Setup_Stream();
        SDR.Start_Stream();

        // debug
        int b = FFTLength/2*ceil(Display.span/Display.resBW);

        for(int z = 0; z<(int)ceil(Display.span/Display.resBW); z++)
        {

            int samplesRead = LMS_RecvStream(&SDR.streamId, StreamBuffer, FFTLength, NULL, 3000);
            usleep(SDRSetting.frameTime*10e6);

            for(int x=0; x<samplesRead; x+=2)
            {
                FFT.in[x/2][0]= StreamBuffer[x];
                FFT.in[x/2][1] = StreamBuffer[x+1];
            }

            FFT.Execute_FFTW();

            for(int y=0; y<FFTLength/2; y++)
            {
                YAxisBuffer[z*FFTLength/2+y] = 10*log10(pow(0.000195313 * pow(FFT.out[y][0]*FFT.out[y][0]+FFT.out[y][1]*FFT.out[y][1], 0.5)/FFTLength, 2)/50/0.001);
                XAxisBuffer[z*FFTLength/2+y] = SDRSetting.centerFreq + (SDRSetting.sampleRate / FFTLength)*y;
            }

            plotFile = fopen("plotFile.txt", "w+");
            for(int x=0; x<(int)ceil(Display.span/Display.resBW)*FFTLength/2; x++)
            {
                fprintf(plotFile, "%f   %f\n", XAxisBuffer[x], YAxisBuffer[x]);
            }
            fclose(plotFile);
            gnuplotHandle.plotfile_xy("plotFile.txt", 1, 2, "Spectrum");

            //clear last plot before making a new one
            gnuplotHandle.reset_plot();

            Change_RF_Settings(SDRSetting.centerFreq+SDRSetting.sampleRate/2 + SDRSetting.sampleRate/4);
        }
    };

    private:

        struct DispSettings
        {
            float_type centerFreq;
            float_type resBW;
            float_type span;
            double refAmp;
        };

        struct SDRSettings
        {
            unsigned gain;
            float_type centerFreq;
            float_type sampleRate;
            float_type frameTime;
        };

        struct PlotSettings
        {
            float_type refAmp;
            FILE* plotFile;
            char* fileName;
        };

        SDRSettings SDRSetting;
        DispSettings Display;
        int16_t *StreamBuffer;
        double *YAxisBuffer;
        double *XAxisBuffer;
        Gnuplot gnuplotHandle = Gnuplot("lines");
        FFT_Settings FFT = FFT_Settings(8192);
        Lime_SDR_mini_Rx SDR = Lime_SDR_mini_Rx();

        int getFFTLength()
        {
            return FFT.DFTlength;
        }

        void Change_RF_Settings(float_type centerFreqIn = 0, float_type sampleRateIn = 0, unsigned gainIn = 100, int FFTIn = 0)
        {
            if(centerFreqIn != 0)
            {
                // change this to respect cases when centerFreq and sampleRate shall be changed
                SDRSetting.centerFreq = centerFreqIn - SDRSetting.sampleRate/4;
                SDR.Change_RF_Settings(SDRSetting.centerFreq);
            }
            if(sampleRateIn != 0)
            {
                SDRSetting.sampleRate = sampleRateIn;
                SDR.Change_RF_Settings(0,SDRSetting.sampleRate);
            }
            if(gainIn != 0)
            {
                SDRSetting.gain = gainIn;
                SDR.Change_RF_Settings(0,0,SDRSetting.gain);
            }
            if(FFTIn != 0)
            {
                FFT.Change_FFT_Settings(FFTIn);
            }
        }
};

int error(char *string, lms_device_t* device)
{
    printf("--------------------------------------------------\n"
           "Error: %s\n"
           "--------------------------------------------------", string);
    if (device != NULL)
        LMS_Close(device);
    exit(-1);
}


//Main Function
int main(int argc, char** argv)
{
    //Initialize variables to be used later on
    float_type CenterFreq = 940e6;
    float_type sample = 10e6;

    Trace Trace1 = Trace(CenterFreq, sample*6, sample, 20, -15, 1024, 1);

    //gnuplotHandle.set_yrange(-150, -5);

    Trace1.Aqcuire_Sweep();


    int x = 0;

    return 0;
}
