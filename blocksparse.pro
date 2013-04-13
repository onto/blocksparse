TEMPLATE = app
CONFIG += console
CONFIG -= qt

SOURCES += \
    main.cpp

HEADERS += \
    sparsematrix.h \
    matrix.h \
    blockmatrix.h \
    matrixoperations.h

QMAKE_CXXFLAGS += -fopenmp
QMAKE_LFLAGS += -fopenmp


