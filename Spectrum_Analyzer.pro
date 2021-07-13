TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += \
        main.c
INCLUDEPATH += /usr/include/c++/9\
               /home/ungare/Documents/Programs/LimeSuite/src/protocols

LIBS += -lLimeSuite -lfftw3

OBJECTS += gnuplot_i.o
