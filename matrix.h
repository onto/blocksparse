#ifndef MATRIX_H
#define MATRIX_H

#include <vector>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <omp.h>

#include "sparsematrix.h"

class Vector
{
public:
    Vector()
    {
        H = 0;
    }
    Vector(const Vector &V1)
    {
        H = V1.H;
        V = V1.V;
    }
    Vector(size_t h)
    {
        H = h;
        V.resize(H+1);
    }

    Vector(std::vector<double> &V1)
    {
        H = V1.size()-1;
        V = V1;
    }

    Vector(const char *file)
    {
        std::ifstream in(file);

        V.push_back(0);

        double x;
        while (in >> x)
            V.push_back(x);

        H = V.size()-1;
    }

    std::vector<double> V;
    size_t H;

    inline double get(size_t i) const
    {
        return V[i];
    }

    inline void set(size_t i, double value)
    {
        V[i] = value;
    }

    void print() const
    {
        for (size_t i = 1; i <= H; ++i)
            std::cout << V[i] << '\n';
    }


    friend Vector operator +(const Vector &V1, const Vector &V2)
    {
        if (V1.H != V2.H) {}

        Vector OV(V1.H);

#pragma omp parallel for
        for (size_t i = 0; i <= V1.H; ++i)
            OV.V[i] = V1.V[i] + V2.V[i];

        return OV;
    }


    friend Vector operator -(const Vector &V1, const Vector &V2)
    {
        if (V1.H != V2.H) {}

        Vector OV(V1.H);

#pragma omp parallel for
        for (size_t i = 0; i <= V1.H; ++i)
            OV.V[i] = V1.V[i] - V2.V[i];

        return OV;
    }


    friend Vector& operator -=(Vector &V1, const Vector &V2)
    {
        if (V1.H != V2.H) {}

#pragma omp parallel for
        for (size_t i = 0; i <= V1.H; ++i)
            V1.V[i] -= V2.V[i];

        return V1;
    }


    friend Vector& operator +=(Vector &V1, const Vector &V2)
    {
        if (V1.H != V2.H) {}

#pragma omp parallel for
        for (size_t i = 0; i <= V1.H; ++i)
            V1.V[i] += V2.V[i];

        return V1;
    }

};

class Matrix
{
public:
    Matrix()
    {
        H = 0;
        W = 0;
    }
    Matrix(const Matrix &_M)
    {
        H = _M.H;
        W = _M.W;
        M = _M.M;
    }
    Matrix(size_t h, size_t w)
    {
        H = h;
        W = w;
        M.resize(H*W+1);
    }
    Matrix(size_t h)
    {
        H = h;
        W = h;
        M.resize(H*W+1);
    }

    std::vector<double> M;
    size_t H;
    size_t W;

    inline double get(size_t row, size_t col) const
    {
        return M[H*(row-1)+col];
    }

    inline void set(size_t row, size_t col, double value)
    {
        M[W*(row-1)+col] = value;
    }

    void print() const
    {
        for (size_t i = 1; i <= H; ++i)
        {
            for (size_t j = 1; j <= W; ++j)
            {
                if (fabs(get(i,j)) < 1e-10)
                    std::cout << "0\t";
                else
                    std::cout << get(i,j) << '\t';
            }
            std::cout << ";\n";
        }
    }

    void swapRow(size_t r1, size_t r2)
    {
        double *t = new double[W];

        memcpy(t, &M[W*(r1-1)+1], W*sizeof(double));
        memcpy(&M[W*(r1-1)+1], &M[W*(r2-1)+1], W*sizeof(double));
        memcpy(&M[W*(r2-1)+1], t, W*sizeof(double));

        delete[] t;
    }

    void swapCol(size_t c1, size_t c2)
    {
#pragma omp parallel for
        for (size_t i = 1; i <= H; ++i)
        {
            size_t h1 = W*(i-1)+c1, h2 = W*(i-1)+c2;
            double t = M[h1];
            M[h1] = M[h2];
            M[h2] = t;
        }
    }


    friend Matrix operator +(const Matrix &M1, const Matrix &M2)
    {
        if ((M1.H != M2.H) || (M1.W != M2.W)) {} // надо бы бросить исключение

        Matrix OM(M1);

#pragma omp parallel for
        for (size_t i = 1; i <= OM.M.size(); ++i)
            OM.M[i] += M2.M[i];

        return OM;
    }


    friend Matrix operator -(const Matrix &M1, const Matrix &M2)
    {
        if ((M1.H != M2.H) || (M1.W != M2.W)) {} // надо бы бросить исключение

        Matrix OM(M1);

#pragma omp parallel for
        for (size_t i = 1; i <= OM.M.size(); ++i)
            OM.M[i] -= M2.M[i];

        return OM;
    }

    friend Matrix& operator -=(Matrix &M1, const Matrix &M2)
    {
        if ((M1.H != M2.H) || (M1.W != M2.W)) {} // надо бы бросить исключение

#pragma omp parallel for
        for (size_t i = 1; i < M1.M.size(); ++i)
            M1.M[i] -= M2.M[i];

        return M1;
    }


    friend Matrix& operator +=(Matrix &M1, const Matrix &M2)
    {
        if ((M1.H != M2.H) || (M1.W != M2.W)) {} // надо бы бросить исключение

#pragma omp parallel for
        for (size_t i = 1; i < M1.M.size(); ++i)
            M1.M[i] += M2.M[i];

        return M1;
    }

    friend Matrix operator *(const Matrix &M1, const Matrix &M2)
    {
        if (M1.W != M2.H) {} // надо бы бросить исключение

        Matrix OM(M1.H, M2.W);

#pragma omp parallel for
        for (size_t i = 1; i <= M1.H; ++i)
        {
            size_t h1 = M1.W*(i-1);
            for (size_t j = 1; j <= M2.W; ++j)
            {
                size_t h2 = OM.W*(i-1)+j;
                for (size_t k = 1; k <= M1.W; ++k)
                    OM.M[h2] += M1.M[h1+k] * M2.get(k,j);
            }
        }

        return OM;
    }


    friend Vector operator *(const Matrix &M, const Vector &V)
    {
        if (M.W != V.H) {} // надо бы бросить исключение

        Vector OV(M.H);

#pragma omp parallel for
        for (size_t i = 1; i <= M.H; ++i)
        {
            size_t h = M.W*(i-1);
            for (size_t j = 1; j <= M.W; ++j)
                OV.V[i] += M.M[h+j] * V.V[j];
        }

        return OV;
    }

};


#endif // MATRIX_H
