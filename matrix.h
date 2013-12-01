#ifndef MATRIX_H
#define MATRIX_H

#include <fstream>
#include <iostream>
#include <vector>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <omp.h>

#include <sparsematrix.h>
//class SparseMatrix;

#define SPARSE_END 0

class Vector
{
public:
    Vector();
    Vector(const Vector &V1);
    Vector(size_t h);
    Vector(std::vector<double> &V1);
    Vector(const char *file);

    std::vector<double> V;
    size_t H;

    inline double get(size_t i) const;
    inline void set(size_t i, double value);

    void zeros(size_t h);
    void resize(size_t h);
    void print() const;
    void permute(std::vector<size_t>& P);

    friend Vector operator +(const Vector &V1, const Vector &V2);
    friend Vector operator -(const Vector &V1, const Vector &V2);

    friend Vector& operator -=(Vector &V1, const Vector &V2);
    friend Vector& operator +=(Vector &V1, const Vector &V2);

    friend Vector operator *(const SparseMatrix &S, const Vector &V);
};

class Matrix
{
public:
    Matrix();
    Matrix(const Matrix &_M);
    Matrix(size_t h, size_t w);
    Matrix(size_t h);
    Matrix(const char *file);
    Matrix(const SparseMatrix &S);

    std::vector<double> M;
    size_t H;
    size_t W;

    inline double get(size_t row, size_t col) const;
    inline void set(size_t row, size_t col, double value);
    void zeros(size_t h, size_t w);
    void print() const;

    void swapRow(size_t r1, size_t r2);
    void swapCol(size_t c1, size_t c2);

    friend Matrix operator +(const Matrix &M1, const Matrix &M2);
    friend Matrix operator -(const Matrix &M1, const Matrix &M2);

    friend Matrix& operator -=(Matrix &M1, const Matrix &M2);
    friend Matrix& operator +=(Matrix &M1, const Matrix &M2);

    friend Matrix operator *(const Matrix &M1, const Matrix &M2);
    friend Vector operator *(const Matrix &M, const Vector &V);

    friend Matrix operator *(const SparseMatrix &S, const Matrix &M);
    friend Matrix operator *(const Matrix &M, const SparseMatrix &S);
    friend Vector operator *(const SparseMatrix &S, const Vector &V);

    friend Matrix operator +(const SparseMatrix &S, const Matrix &M);
    friend Matrix operator +(const Matrix &M, const SparseMatrix &S);

    friend Matrix& operator +=(Matrix &M, const SparseMatrix &S);

    friend Matrix operator -(const SparseMatrix &S, const Matrix &M);
    friend Matrix operator -(const Matrix &M, const SparseMatrix &S);
};

#endif // MATRIX_H
