TEMPLATE = app
CONFIG += console
CONFIG -= qt

SOURCES += generator.cpp

QMAKE_CXXFLAGS += -fopenmp
QMAKE_LFLAGS += -fopenmp

