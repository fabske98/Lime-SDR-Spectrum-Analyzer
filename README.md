# Lime-SDR-Spectrum-Analyzer

A user configurable Spectrum Analyzer Application featuring a frequency range from 10 MHz to 3.5 GHz and a maximum resolution bandwidth of 15 MHz using a Lime SDR (Software-Defined-Radio) Mini. QT Creator IDE has been used.

See an example spectrum Center Frequency 960 MHz, Span 40MHz, Resolution Bandwidth 10 MHz, Clear Write 
![Spectrum](https://github.com/fabske98/Lime-SDR-Spectrum-Analyzer/blob/C%2B%2B_OOP_Modular/images/Snapshot_No_Averaging.png)

See an example spectrum Center Frequency 960 MHz, Span 40MHz, Resolution Bandwidth 10 MHz, Averaging Factor 30
![Spectrum](https://github.com/fabske98/Lime-SDR-Spectrum-Analyzer/blob/C%2B%2B_OOP_Modular/images/Snapshot_Averaging_30.png)

See the hardware setup
![Hardware](https://github.com/fabske98/Lime-SDR-Spectrum-Analyzer/blob/C%2B%2B_OOP_Modular/images/Hardware.jpg)

This project is still under development. Current features include

- Setting the following common Spectrum Analyzer settings upon creating a trace Center Frequency, Span, Resolution Bandwidth, Gain, Reference Amplitude and Trace averaging.
- Plotting the Spectrum (using Gnuplot)
- Creating an own trace for each LimeSDRmini connected

The following topics are to be tackled soon

- Remove the DC bias and center spike (caused by the SDR) from the spectrum
- Changing of Spectrum Analyzer settings on the fly (right now the created Trace object needs to be destroyed and a new to be created)
- Use the (internal) second Rx channel for frequencies greater than 2 GHz
- User defineable and automatically (e.g. max) set markers in the plot
- Channel power measurements
- A GUI with a similar user interface as one is accustomed to from spectrum analyzers
- Release as ver. 1.0 standalone application once the topics above are implemented

## Installation

The following external libraries  sources are used for this project and are required for using the code published in this repository

- limesuite (SDR driver) [Link](httpsgithub.commyriadrfLimeSuite)
- gnuplot_i (Plotting) (sourcefile and header already include in this repo)
- fftw3 (Fast-Fourier-Transformation) [Link](httpfftw.org)

```bash
pip install foobar
```

## Usage

```cpp
#include fftw3.h
#include limesuitelimesuite.h
#include gnuplot_i.h

Trace myTrace();
myTrace.Aqcuire_Sweep();

```

## Contributing
Pull requests are welcome. For major changes, please contact me under fabian.toeroek98@gmail.com to discuss your ideas.

## License
[GPLv3](httpschoosealicense.comlicensesgpl-3.0)
