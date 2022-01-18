/*
Description:
A simple spectrum analyzer application using a LimeSDR mini, FFTW and GNU Plot.

Author:
Fabian Török

V1.0
 */
//**************************************************************************************
//#Project Includes
//**************************************************************************************
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include "lime/LimeSuite.h"     // Required to drive LimeSDR
#include <gnuplot_i.hpp>          // Required for Plotting
#include <fftw3.h>              // Required for DFT and FFT
#include <unistd.h>

//**************************************************************************************
//#Function:    int error(char*, lms_device_t*)
//**************************************************************************************
//#Description: This function outputs an error message and closes the connection to the
//              SDR.
//**************************************************************************************
//#Parameters:  <string>    :[---] <char*> Error Message to be displayed
//              <device>    :[---] <lms_device_t*> Device handle
//**************************************************************************************
//#Date         18.01.2022
//**************************************************************************************
//#Author       Fabian Török
//**************************************************************************************
int error(char *string, lms_device_t* device)
{
    printf("--------------------------------------------------\n"
           "Error: %s\n"
           "--------------------------------------------------", string);
    if (device != NULL)
        LMS_Close(device);
    exit(-1);
}

//**************************************************************************************
//#Class:       Lime_SDR_mini_Rx
//**************************************************************************************
//#Description: This class manages the driver interaction with the SDR. All settings are
//              set through the Setup_SDR or Change_SDR_Settings functions.
//**************************************************************************************
//#Members:     Public
//#Variables
//              <device>                :[---] <lms_device_t*> Device Handle
//              <streamId>              :[---] <lms_stream_t> Stream handle
//#Functions
//              <Lime_SDR_mini_Rx>      :Constructor
//              <~Lime_SDR_mini_Rx>     :Destructor
//              <Setup_SDR>             :Opens SDR connection and sets specified
//                                       settings
//              <Change_RF_Settings>    :Changes RF Settings
//              <Setup_Stream>          :Stream Setup
//              <Start_Stream>          :Start Stream
//              <Stop_Stream>           :Stop Stream
//              <Close_Stream>          :Destroy Stream
//**************************************************************************************
//#Date         18.01.2022
//**************************************************************************************
//#Author       Fabian Török
//**************************************************************************************
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

            Setup_Stream();
            Start_Stream();

            return 1;
        }

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
        }

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
        }

        void Start_Stream()
        {
            //Start streaming
            LMS_StartStream(&streamId);
        }

        void Stop_Stream()
        {
            //Stop streaming
            LMS_StopStream(&streamId);
        }

        void Close_Stream()
        {
            LMS_DestroyStream(device, &streamId);
        }
};

//**************************************************************************************
//#Class:       FFT_Settings
//**************************************************************************************
//#Description: This class manages the Fast Fourier Transformations. The FFT Length is
//              only changeable through this class
//**************************************************************************************
//#Members:     Public
//#Variables
//              <DFTlength>             :[---] <int> Fourier Transform Length
//              <in>                    :[---] <fftw_complex*>(double*) Input to FFT
//              <out>                   :[---] <fftw_complex*>(double*) Output from FFT
//#Functions
//              <FFT_Settings>          :Constructor
//              <~FFT_Settings>         :Destructor
//              <Execute_FFTW>          :Performs FFT with <in> and saves result to <out>
//              <Change_FFT_Settings>   :Changes FFT Length and adapts size of <in> and
//                                       <out>
//              Private
//#Variables
//              <p>                     :<fftw_plan> this structure contains FFT settings
//**************************************************************************************
//#Date         18.01.2022
//**************************************************************************************
//#Author       Fabian Török
//**************************************************************************************
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
            // allocate fftw_complex array of the size DFTlength * 2
            in = (fftw_complex*) fftw_malloc(DFTlength * sizeof(fftw_complex));
            out = (fftw_complex*) fftw_malloc(DFTlength * sizeof(fftw_complex));
            // declare FFTLength, input and output as well as Forward FFT
            p = fftw_plan_dft_1d(DFTlength, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
        };

        void Execute_FFTW()
        {
            fftw_execute(p);
        };

        void Change_FFT_Settings(int DFTlengthIn)
        {
            // Free old memory
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
            // Free dynamically allocated memory
            fftw_destroy_plan(p);
            fftw_free(in);
            fftw_free(out);
        };
};

//**************************************************************************************
//#Class:       FFT_Settings
//**************************************************************************************
//#Description: This class manages the Fast Fourier Transformations. The FFT Length is
//              only changeable through this class
//**************************************************************************************
//#Members:     Public
//#Functions
//              <Trace>                 :Constructor
//              <~Trace>                :Destructor
//              <Aqcuire_Sweep>         :Performs a sweep. If user defined span is
//                                       greater than the resolution bandwidth, multiple
//                                       frames of <FFTLength> are aquired after retuning
//                                       the SDR
//              <Plot_Spectrum>         :Plots content of XAxisBuffer and YAxisBuffer
//
//              Private
//#Variables
//              <SDRSetting>            :<SDRSettings> struct of SDR settings
//              <DispSetting>           :<DisplaySettings> struct of RF settings of
//                                       Spectrum Analyzer
//              <PlotSetting>           :<PlotSettings> struct of Plot specific settings
//              <StreamBuffer>          :<int16_t*> array for storing stream data
//              <XAxisbuffer>           :<double*> array for storing XAxis plot data
//              <YAxisbuffer>           :<double*> array for storing YAxis plot data
//              <YAxisbufferAvg>        :<double*> array for storing YAxis plot for
//                                       averaging
//              <GnuplotHandle>         :<Gnuplot> Gnuplot object
//              <FFT>                   :<FFT_Settings> FFT_Settings object
//              <SDR>                   :<Lime_SDR_mini_Rx> Lime_SDR_mini_Rx object
//              <AveragingCounter>      :<int> counter to track saved frames for
//                                       averaging
//              <AveragingCounter>      :<bool> flag indicating that first averaging
//                                       cycle has finished
//#Functions
//              <getFFTLength>          :returns current FFT length
//              <Change_RF_Settings>    :translates user requirements and changes SDR
//                                       RF settings as well
//
//**************************************************************************************
//#Date         18.01.2022
//**************************************************************************************
//#Author       Fabian Török
//**************************************************************************************
class Trace
{
    public:
        Trace(float_type centerFreqIn, float_type dispSpanIn,  float_type resBWIn, double gainIn, double refAmpIn, int FFTIn, int averaging, unsigned TraceNum):FFT(8192)
        {
            // Apply user settings
            Display.centerFreq = centerFreqIn;
            Display.span = dispSpanIn;
            Display.resBW = resBWIn;
            Display.refAmp = refAmpIn;
            Display.averaging = averaging;
            initAveraging = false;
            AveragingCounter = 1;
            // Specify SDR sample rate according to Nyquist Theorem
            SDRSetting.sampleRate = resBWIn*2;
            // Calculate time for SDR to sample a frame
            SDRSetting.frameTime = 1/SDRSetting.sampleRate*FFTIn;
            // Set initial Center Frequency to beginning of span. Note: the term Center Frequency is misleading here as the frequency set as "Center Frequency"
            // in the SDR will from where the downmixed Basband signal starts. ->Start Frequency is a better name here
            SDRSetting.centerFreq = centerFreqIn - dispSpanIn / 2;
            SDRSetting.gain = gainIn;
            FFT.Change_FFT_Settings(FFTIn);
            // Length of the dynaminc arrays equals FFTIn/2*(Display.span/Display.resBW --> FFTIn/2 is sufficient here as the other half of the FFT will be redundant
            // information. See Nyquist Theorem.
            YAxisBuffer = new double [int(ceil(FFTIn/2*(Display.span/Display.resBW)))];
            YAxisBufferAvg = new double [int(ceil(FFTIn/2*(Display.span/Display.resBW)))*averaging];
            XAxisBuffer = new double [int(ceil(FFTIn/2*(Display.span/Display.resBW)))];
            for(int x=0; x<(int)(ceil(FFTIn/2*(Display.span/Display.resBW))); x++)
                YAxisBuffer[x] = 0;
            // The averaging buffer will be <averaging> times bigger than the rest
            for(int x=0; x<(int)(ceil(FFTIn/2*(Display.span/Display.resBW))*averaging); x++)
                YAxisBufferAvg[x] = 0;
            for(int x=0; x<int(ceil(FFTIn/2*(Display.span/Display.resBW))); x++)
                XAxisBuffer[x] = SDRSetting.centerFreq + (SDRSetting.sampleRate / FFT.DFTlength)*x;
            // Stream will fit the data as follows I,Q,I,Q,I,Q... --> double the size of DFTlength
            StreamBuffer = new int16_t[FFT.DFTlength * 2];
            // Setup SDR
            SDR.Setup_SDR(SDRSetting.centerFreq, SDRSetting.sampleRate, SDRSetting.gain);
        };
        //**************************************************************************************
        //**************************************************************************************
        ~Trace()
        {
            // free dynamically allocated arrays
            free(YAxisBuffer);
            free(XAxisBuffer);
            free(StreamBuffer);
            free(YAxisBufferAvg);
        };

        void Aqcuire_Sweep()
        {
            // get current FFTLength
            int FFTLength = getFFTLength();
            // calculate amount of bins in one sweep (from StartFreq until StopFreq)
            int sweepBins = (int)ceil(Display.span/Display.resBW*FFTLength/2);

            // repeat this for Span/resBW times
            for(int z = 0; z<(int)ceil(Display.span/Display.resBW); z++)
            {
                // receive samples from SDR
                int samplesRead = LMS_RecvStream(&SDR.streamId, StreamBuffer, FFTLength, NULL, 3000);
                // sleep the time required for samples to be received
                usleep(SDRSetting.frameTime*10e6);

                // Fill FFT.in --> I data goes into FFT.in[x][0], Q data into FFT.in[x][1]
                for(int x=0; x<samplesRead; x+=2)
                {
                    FFT.in[x/2][0]= StreamBuffer[x];
                    FFT.in[x/2][1] = StreamBuffer[x+1];
                }

                // Calculate FFT
                FFT.Execute_FFTW();

                // If averaging is not required
                if(Display.averaging == 1)
                {
                    for(int y=0; y<FFTLength/2; y++)
                    {   // 0.000195313 V equals one step of ADC (only valid when stream data is configured as 12-Bit integers!!)
                        // Calulate Amplitude of signal by (I²+Q²)¹/² = A, (divide by FFTLength!) --> then calculate Power P = A²/R --> then convert to logarithmic scale
                        YAxisBuffer[z*FFTLength/2+y] = 10*log10(pow(0.000195313 * pow(FFT.out[y][0]*FFT.out[y][0]+FFT.out[y][1]*FFT.out[y][1], 0.5)/FFTLength, 2)/50/0.001);
                        //XAxisBuffer[z*FFTLength/2+y] = SDRSetting.centerFreq + (SDRSetting.sampleRate / FFTLength)*y;
                    }
                }
                // Averaging required
                else
                {
                    for(int y=0; y<FFTLength/2; y++)
                    {
                        // Paste newly aqcuired samples into the averaging buffer
                        YAxisBufferAvg[AveragingCounter*sweepBins+(z*FFTLength/2+y)] = pow(0.000195313 * pow(FFT.out[y][0]*FFT.out[y][0]+FFT.out[y][1]*FFT.out[y][1], 0.5)/FFTLength, 2)/50;

                        if(initAveraging == false)
                        {
                            double averagedPower = 0;
                            // Add up power for this bin from current + previous sweeps
                            for(int x=0; x<AveragingCounter+1; x++)
                                averagedPower += YAxisBufferAvg[sweepBins*x+(z*FFTLength/2+y)];

                            // Divide by number of previous + current sweeps stored
                            averagedPower /= AveragingCounter+1;

                            // convert to dBm
                            YAxisBuffer[z*FFTLength/2+y] = 10*log10(averagedPower/0.001);

                            // set initaveraging to true, once enough sweeps are saved for averaging as specified by user
                            if(AveragingCounter == Display.averaging-1)
                                initAveraging = true;
                        }
                        else
                        {
                            double averagedPower = 0;
                            // Add up power for this bin from current + previous sweeps
                            for(int x=0; x<Display.averaging; x++)
                                averagedPower += YAxisBufferAvg[sweepBins*x+(z*FFTLength/2+y)];
                            // Divide by number of previous + current sweeps stored
                            averagedPower /= Display.averaging;

                            // Convert to logarithmic scale and assign to YAxisBuffer
                            YAxisBuffer[z*FFTLength/2+y] = 10*log10(averagedPower/0.001);
                        }
                    }
                }
                // Plot the spectrum
                Plot_Spectrum();

                // Tune SDR for next step in the sweep
                Change_RF_Settings(SDRSetting.centerFreq+SDRSetting.sampleRate/2 + SDRSetting.sampleRate/4);
            }

            if(AveragingCounter == Display.averaging-1)
                AveragingCounter = 0;
            else
                AveragingCounter += 1;
            Change_RF_Settings(Display.centerFreq-Display.span/2+SDRSetting.sampleRate/4);
        };
        //**************************************************************************************
        //**************************************************************************************
        void Plot_Spectrum()
        {
            FILE* plotFile;
            int FFTLength = getFFTLength();
            // set upper limit of plot to refAmp
            GnuplotHandle.set_yrange(-150, Display.refAmp);

            // write plotFile from contents of XAxisBuffer and YAxisBuffer
            plotFile = fopen("plotFile.txt", "w+");
            for(int x=0; x<(int)ceil(Display.span/Display.resBW)*FFTLength/2; x++)
            {
                fprintf(plotFile, "%f   %f\n", XAxisBuffer[x], YAxisBuffer[x]);
            }
            fclose(plotFile);

            GnuplotHandle.plotfile_xy("plotFile.txt", 1, 2, "Spectrum");

            // set bounding box for channel measurement
            //GnuplotHandle.cmd("set object rectangle from screen 0.4,0.06 to screen 0.5,0.97 behind fillcolor rgb 'red' fillstyle solid 0.5 border");

            //clear last plot before making a new one
            GnuplotHandle.reset_plot();
        };
        //**************************************************************************************
        //**************************************************************************************

    private:
        // Diplay settings struct
        struct DispSettings
        {
            float_type centerFreq;
            float_type resBW;
            float_type span;
            double refAmp;
            int averaging;
        };
        // SDR settings struct
        struct SDRSettings
        {
            unsigned gain;
            float_type centerFreq;
            float_type sampleRate;
            float_type frameTime;
        };
        // Plot settings struct
        struct PlotSettings
        {
            float_type refAmp;
            FILE* plotFile;
            char* fileName;
        };

        // declare private variables
        SDRSettings SDRSetting;
        DispSettings Display;
        int16_t *StreamBuffer;
        double *YAxisBuffer;
        double *YAxisBufferAvg;
        double *XAxisBuffer;
        Gnuplot GnuplotHandle = Gnuplot("lines");
        FFT_Settings FFT = FFT_Settings(8192);
        Lime_SDR_mini_Rx SDR = Lime_SDR_mini_Rx();
        int AveragingCounter;
        bool initAveraging;

        // Returns FFT Length
        int getFFTLength()
        {
            return FFT.DFTlength;
        };

        // Adapt user input to correct settings for SDR
        void Change_RF_Settings(float_type centerFreqIn=0, float_type sampleRateIn=0, unsigned gainIn=0, int FFTIn=0)
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
        };
};


//Main Function
int main(int argc, char** argv)
{
    //Initialize variables to be used later on
    float_type centerFreq = 965e6;
    float_type span = 40e6;
    float_type sample = 10e6;
    double gain = 25;
    double refAmp = -15;
    int FFTlength = 1024;
    int averaging = 30;
    unsigned nTrace = 1;

    Trace Trace1 = Trace(centerFreq, span, sample, gain, refAmp, FFTlength, averaging, nTrace);

    // Aqcuire 100 sweeps
    for(int u=0; u<100; u++)
        Trace1.Aqcuire_Sweep();

    return 0;
}

