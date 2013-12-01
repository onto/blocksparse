TEMPLATE = app
CONFIG += console
CONFIG -= qt

SOURCES += \
    main.cpp \
    matrix.cpp \
    sparsematrix.cpp

HEADERS += \
    sparsematrix.h \
    matrix.h \
    blockmatrix.h \
    matrixoperations.h \
    matrixcontainers.h \
    cudaoperations.h \
    decompositor.h \
    udecompositor.h

QMAKE_CXXFLAGS += -Wextra -fopenmp
QMAKE_LFLAGS += -fopenmp

