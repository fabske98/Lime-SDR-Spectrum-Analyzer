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
#include <gnuplot_i.h>          // Required for Plotting
#include <time.h>
#include <math.h>
#include <fftw3.h>              // Required for DFT and FFT
#include <csvparser.h>
#ifdef USE_GNU_PLOT
#include "gnuPlotPipe.h"
#endif
#define CONFIG_UPDATE_TIME 1000



//Device structure, should be initialized to NULL
lms_device_t* device = NULL;

int error(char *string)
{
    printf("--------------------------------------------------\n"
           "Error: %s\n"
           "--------------------------------------------------", string);
    if (device != NULL)
        LMS_Close(device);
    exit(-1);
}


class Lime_SDR_mini_Rx
{
    public:
        lms_device_t* device;
        lms_stream_t streamId;

    private:
        Lime_SDR_mini_Rx()
        {
            device = NULL;
        };
        template <class Disp_Meas_Settings>
        void setup_SDR(Disp_Meas_Settings Settings_In)
        {
            //open the first device
            lms_info_str_t list[8];
            if (LMS_Open(&device, list[0], NULL))
                error("Alio");

            if (LMS_Init(device) != 0)
                error("Alio");

            //Enable RX channel
            //Channels are numbered starting at 0
            if (LMS_EnableChannel(device, LMS_CH_RX, 0, true) != 0)
                error("Alio");

            //Set center frequency
            if (LMS_SetLOFrequency(device, LMS_CH_RX, 0, Settings_In.Center_Freq) != 0)
                error("Alio");

            //Set sample rate, ask to use 2x oversampling in RF
            //This sets sampling rate for all channels
            if (LMS_SetSampleRate(device, Settings_In.Sample_Rate, 1) != 0)
                error("Alio");

            //Set Lowpass Filter Bandwidth
            LMS_SetLPFBW(device, LMS_CH_RX, 0, Settings_In.Sample_Rate*2);

            LMS_SetGaindB(device, LMS_CH_RX, 0, Settings_In.Gain);
        };
        void Setup_Stream()
        {
            streamId.channel = 0;               //channel number
            streamId.fifoSize = 1024 * 1024;    //fifo size in samples
            streamId.throughputVsLatency = 1.0; //optimize for max throughput
            streamId.isTx = false;              //RX channel
            streamId.dataFmt = lms_stream_t::LMS_FMT_I12;     //12-bit integers //lms_stream_t::

            //Setup Stream
            if (LMS_SetupStream(device, &streamId) != 0)
                error("Alio");
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
        ~Lime_SDR_mini_Rx()
        {
            LMS_Close(device);
        };
};

class Disp_Meas_Settings
{
    protected:
        float_type Disp_CenterFreq;
        float_type Disp_ResBw;
        float_type Sweep_Time;
        double Gain;
        double LPBWF_BW;
        double Ref_Amp;
        float_type Center_Freq;
        float_type Sample_Rate;

    public:
        Disp_Meas_Settings(float_type Center_Freq_In, float_type Sample_Rate_In, double Gain_In, double LPBWF_BW_In, double Ref_Amp_In)
        {
            Center_Freq = Center_Freq_In - Sample_Rate_In/4;
            Sample_Rate = Sample_Rate_In;
            Gain = Gain_In;
            LPBWF_BW = LPBWF_BW_In;
            Ref_Amp = Ref_Amp_In;
        }

        void Change_RF_Settings(float_type Center_Freq_In, float_type Sample_Rate_In)
        {
            Center_Freq = Center_Freq_In - Sample_Rate_In/4;
            Sample_Rate = Sample_Rate_In;
        }

        Disp_Meas_Settings();
};

class FFT_Settings
{
    protected:
        int DFTlength;
        double *PSD_1Hz;
        double *FREQ;
        double *tempArr;
        fftw_complex *in;
        fftw_complex *out;
        fftw_plan p;

    public:
        FFT_Settings(int DFTlengthIn, double gainIn)
        {
            DFTlength = DFTlengthIn;
            PSD_1Hz = new double [DFTlength/2];
            FREQ = new double [DFTlength/2];
            in = new fftw_complex [DFTlength];
            out = new fftw_complex [DFTlength];
            p = fftw_plan_dft_1d(DFTlength, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
        };
        void Execute_FFTW()
        {
            fftw_execute(p);
        };
        void Change_FFT_Settings(int DFTlengthIn)
        {
            //Free old memory
            free(PSD_1Hz);
            free(FREQ);
            fftw_free(in);
            fftw_free(out);
            fftw_destroy_plan(p);

            // Set new values
            DFTlength = DFTlengthIn;
            PSD_1Hz = new double [DFTlength/2];
            FREQ = new double [DFTlength/2];
            in = new fftw_complex [DFTlength];
            out = new fftw_complex [DFTlength];
            p = fftw_plan_dft_1d(DFTlength, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
        };


        ~FFT_Settings()
        {
            free(PSD_1Hz);
            free(FREQ);
            fftw_free(in);
            fftw_free(out);
            fftw_destroy_plan(p);
        };
};

/*
int loadConfig()
{
    FILE *fid;
    fid = fopen("config.csv", "r");

    if(fid == NULL)
    {
        error("config file not found");
    }
    else
    {
        fclose(fid);
        CsvParser *csvparser = CsvParser_new("config.csv", ",", 0);
        CsvRow *row;
        while ((row = CsvParser_getRow(csvparser)) )
        {
            printf("NEW LINE:\n");
            const char **rowFields = CsvParser_getFields(row);
            for (int i = 0 ; i < CsvParser_getNumFields(row) ; i++) {
                printf("FIELD: %s\n", rowFields[i]);
            }
            CsvParser_destroy_row(row);
        }
        CsvParser_destroy(csvparser);
        return 0;
    }
}
*/


//Main Function
int main(int argc, char** argv)
{
    //Initialize variables to be used later on
    gnuplot_ctrl *gnuplotHandle;
    int DFTlength = 8192;           //DFT Length
    float_type centerFreq = 433.9e6;  //Center Frequency
    fftw_complex in[DFTlength];     //Input to fftw
    fftw_complex out[DFTlength];    //Output of fftw
    float_type sampleRate = 1e6;    //Bandwidth
    centerFreq = centerFreq - sampleRate/4;
    double PSD_1Hz[DFTlength/2];          //Array for saving power values
    double FREQ[DFTlength/2];         //Array for frequency power values
    double tempArr[DFTlength/2];

    //Setup Stream
    if (LMS_SetupStream(device, &streamId) != 0)
        error("Alio");

    //Initialize data buffers
    const int sampleCnt = DFTlength; //complex samples per buffer
    int16_t buffer[sampleCnt * 2]; //buffer to hold complex values (2*samples))

    //Start streaming
    LMS_StartStream(&streamId);

    //Set this variable to 1 to break out of while loop
    int breakLoop=0;

    //Initialize gnuplot
    gnuplotHandle = gnuplot_init();
    gnuplot_setstyle(gnuplotHandle, "lines");
    gnuplot_cmd(gnuplotHandle, "set yrange [%i:%f]", -150, -30);

    //Get current chip temperature
    LMS_GetChipTemperature(device, 0, tempe);

    //Repeat loop as long as chip temperature remains lower than 55°C and break condition is set to 0
    while(*tempe < 55 && breakLoop==0)
    {
        //Get current chip temperature
        LMS_GetChipTemperature(device, 0, tempe);

        //Receive samples and store number of samples
        int samplesRead = LMS_RecvStream(&streamId, buffer, sampleCnt, NULL, 1000);
        //I and Q samples are interleaved in buffer: IQIQIQ...

        //Store In-Phase values into in[x][0]
        //Store Quadrature values into in[x][1]
        for(int x=0; x<sampleCnt*2; x+=2)
        {
            in[x/2][0]= buffer[x];
            in[x/2][1] = buffer[x+1];
        }

        //Execute FFT
        fftw_execute(p);

        // Calculate bin width and start frequency
        double binWidth = sampleRate/DFTlength;
        double binStart = centerFreq;

        //calculate Power per Hz from FFT --> First calculate Amplitude, then Power (50 Ohm System)
        for (int y=0; y<DFTlength/2; y++)
        {
            PSD_1Hz[y] = 10*log10(pow(0.000195313 * pow(out[y][0]*out[y][0]+out[y][1]*out[y][1], 0.5)/DFTlength*2, 2)/50/0.001); // Try to resemble x * 0.000000429153V/Hz here
            tempArr[y] = PSD_1Hz[y];
            // Shift Array
        }

        //calculate frequency positions after DFT
        for(int y = 0; y<DFTlength/2; y++)
        {
            FREQ[y] = binStart + y*binWidth;
        }

        //plot spectrum
        gnuplot_plot_xy(gnuplotHandle, FREQ, PSD_1Hz, DFTlength/2, "Spec");

        //delay in order to avoid overload of gnuplot
        for(int c = 0; c<5000; c++)
        {
            for(int d = 0; d<5000; d++)
            {
            }
        }

        //clear last plot before making a new one
        gnuplot_resetplot(gnuplotHandle);
    }


    //Stop streaming
    LMS_StopStream(&streamId); //stream is stopped but can be started again with LMS_StartStream()
    LMS_DestroyStream(device, &streamId); //stream is deallocated and can no longer be used

    //Close device
    LMS_Close(device);

    //free allocated memory
    free(gnuplotHandle);
    free(tempe);
    free(gainPtr);
    fftw_destroy_plan(p);


    return 0;
}
