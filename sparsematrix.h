#ifndef SPARSEMATRIX_H
#define SPARSEMATRIX_H

#include <vector>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <omp.h>

//#include "matrix.h"

#define SPARSE_END 0

class SparseMatrix
{
public:
    SparseMatrix();
    SparseMatrix(size_t h, size_t w);
    SparseMatrix(const SparseMatrix &S);
    SparseMatrix(const char *file);

    std::vector<double> V;  //values
    std::vector<size_t> N;  //next in row
    std::vector<size_t> C;  //column
    std::vector<size_t> F;  //first in row

    size_t W; //width
    size_t H; //height

    void set(size_t row, size_t col, double value);
    double get(size_t row, size_t col) const;
    void rm(size_t row, size_t col);
    void add(size_t row, size_t col, double value);

    void swapCol(size_t c1, size_t c2);
    void swapRow(size_t r1, size_t r2);

    void permute(std::vector<size_t>& P);
    void permute(std::vector<size_t>& Pr, std::vector<size_t>& Pc);

    void print() const;
    void save2file(const char *file) const;
    void save2fileold(const char *file) const;
};

#endif // SPARSEMATRIX_H
