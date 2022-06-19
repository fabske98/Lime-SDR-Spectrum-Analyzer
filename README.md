# Lime-SDR-Spectrum-Analyzer

A user configurable Spectrum Analyzer Application featuring a frequency range from 10 MHz to 3.5 GHz and a maximum resolution bandwidth of 15 MHz using a Lime SDR (Software-Defined-Radio) Mini. QT Creator IDE has been used.

See an example spectrum Center Frequency: 960 MHz, Span: 40MHz, Resolution Bandwidth: 10 MHz, Clear Write 
![Spectrum](https://github.com/fabske98/Lime-SDR-Spectrum-Analyzer/blob/Measurement_Options/images/Snapshot_No_Averaging.png)

See an example spectrum Center Frequency: 960 MHz, Span: 40MHz, Resolution Bandwidth: 10 MHz, Averaging Factor 30, Channel Bandwidth Measurement with 20 MHz BW and Marker placed
![Spectrum](https://github.com/fabske98/Lime-SDR-Spectrum-Analyzer/blob/Measurement_Options/images/Snapshot_Averaging_30.png)

See the hardware setup
![Hardware](https://github.com/fabske98/Lime-SDR-Spectrum-Analyzer/blob/C%2B%2B_OOP_Modular/images/Hardware.jpg)

This project is still under development. Current features include

- Setting the following common Spectrum Analyzer settings upon creating a trace Center Frequency, Span, Resolution Bandwidth, Gain, Reference Amplitude and Trace averaging.
- Plotting the Spectrum (using Gnuplot)
- Creating an own trace for each LimeSDRmini connected
- Setting up to 3 markers in the plot
- Performing a Channel Power Measurement

The following topics are to be tackled soon

- Remove the DC bias and center spike (caused by the SDR) from the spectrum
- Changing of Spectrum Analyzer settings on the fly (right now the created Trace object needs to be destroyed and a new to be created)
- Use the (internal) second Rx channel for frequencies greater than 2 GHz
- A GUI with a similar user interface as one is accustomed to from spectrum analyzers
- Release as ver. 1.0 standalone application once the topics above are implemented

## Installation

The following external libraries  sources are used for this project and are required for using the code published in this repository

- limesuite (SDR driver) [Link](https://github.com/myriadrf/LimeSuite)
- gnuplot_i (Plotting API) (sourcefile and header already include in this repo)
- Gnuplot [Link](https://sourceforge.net/projects/gnuplot/files/gnuplot/) (Version 5.4.3 is tested and working)
- fftw3 (Fast-Fourier-Transformation) [Link](http://fftw.org)

- Run and modify the code in Lime_SDR_Spectrum_Analyzer.cpp to your liking
## Usage

```cpp
#include fftw3.h
#include limesuitelimesuite.h
#include gnuplot_i.h

float_type centerFreq = 965e6;
float_type span = 40e6;
float_type sample = 10e6;
double gain = 25;
double refAmp = -15;
int FFTlength = 1024;
int averaging = 30;
unsigned nTrace = 1;

Trace Trace1 = Trace(centerFreq, span, sample, gain, refAmp, FFTlength, averaging, nTrace);

Trace.Aqcuire_Sweep();

```

## Contributing
Pull requests are welcome. For major changes, please contact me under fabian.toeroek98@gmail.com to discuss your ideas.

## License
[GPLv3](httpschoosealicense.comlicensesgpl-3.0)
