TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += \
    Lime_SDR_Spectrum_Analyzer.cpp
INCLUDEPATH += /usr/include/c++/9\
               /home/ungare/Documents/Programs/LimeSuite/src/protocols

LIBS += -lLimeSuite -lfftw3

OBJECTS += gnuplot_i.o

HEADERS += \
    gnuplot_i.hpp
