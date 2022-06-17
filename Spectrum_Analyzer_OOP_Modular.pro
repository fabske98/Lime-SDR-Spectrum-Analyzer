TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += \
    Lime_SDR_Spectrum_Analyzer.cpp

INCLUDEPATH += /usr/include/c++/9\
               /home/ungare/Documents/libraries/LimeSuite/src
INCLUDEPATH += /usr/include/c++/9\
               /home/ungare/Documents/libraries/fftw-3.3.10/api

LIBS += -lLimeSuite -lfftw3

OBJECTS += gnuplot_i.o

HEADERS += \
    gnuplot_i.hpp

win32:CONFIG(release, debug|release): LIBS += -L$$PWD/../../libraries/LimeSuite/build/src/release/ -lLimeSuite
else:win32:CONFIG(debug, debug|release): LIBS += -L$$PWD/../../libraries/LimeSuite/build/src/debug/ -lLimeSuite
else:unix: LIBS += -L$$PWD/../../libraries/LimeSuite/build/src/ -lLimeSuite

INCLUDEPATH += $$PWD/../../libraries/LimeSuite/build/src
DEPENDPATH += $$PWD/../../libraries/LimeSuite/build/src

DISTFILES += \
    gnuplot_i.o
