TEMPLATE = app
CONFIG += console
CONFIG -= qt

SOURCES += \
    main.cpp

HEADERS += \
    sparsematrix.h \
    matrix.h \
    blockmatrix.h \
    matrixoperations.h \
    matrixcontainers.h \
    cudaoperations.h \
    decompositor.h

QMAKE_CXXFLAGS += -fopenmp -Wextra
QMAKE_LFLAGS += -fopenmp

