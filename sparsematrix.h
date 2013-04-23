#ifndef SPARSEMATRIX_H
#define SPARSEMATRIX_H

#include <vector>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <omp.h>

#include "matrix.h"

#define SPARSE_END 0

class SparseMatrix
{
public:
    SparseMatrix();
    SparseMatrix(size_t h, size_t w);
    SparseMatrix(const SparseMatrix &S);
    SparseMatrix(const char *file);

    std::vector<double> V;       //values
    std::vector<size_t> N;  //next in row
    std::vector<size_t> C;  //column
    std::vector<size_t> F;  //first in row

    size_t W; //width
    size_t H; //height

    void set(size_t row, size_t col, double value);
    double    get(size_t row, size_t col) const;
    void rm(size_t row, size_t col);
    void add(size_t row, size_t col, double value);

    void swapCol(size_t c1, size_t c2);
    void swapRow(size_t r1, size_t r2);

    void print() const;
    void save2file(const char *file) const;
    void save2fileold(const char *file) const;


    friend Matrix operator *(const SparseMatrix &S, const Matrix &M);
    friend Matrix operator *(const Matrix &M, const SparseMatrix &S);
    friend Vector operator *(const SparseMatrix &S, const Vector &V);

    friend Matrix operator +(const SparseMatrix &S, const Matrix &M);
    friend Matrix operator +(const Matrix &M, const SparseMatrix &S);

    friend Matrix& operator +=(Matrix &M, const SparseMatrix &S);

    friend Matrix operator -(const SparseMatrix &S, const Matrix &M);
    friend Matrix operator -(const Matrix &M, const SparseMatrix &S);

};



SparseMatrix::SparseMatrix()
{
    H = 0;
    W = 0;

    F.push_back(SPARSE_END);
    V.push_back(0);
    C.push_back(0);
    N.push_back(SPARSE_END);
}


SparseMatrix::SparseMatrix(size_t h, size_t w)
{
    H = h;
    W = w;

    F.assign(h+1, SPARSE_END);
    V.push_back(0);
    C.push_back(0);
    N.push_back(SPARSE_END);
}


SparseMatrix::SparseMatrix(const SparseMatrix &S)
{
    V = S.V;
    N = S.N;
    C = S.C;
    F = S.F;
    H = S.H;
    W = S.W;
}


SparseMatrix::SparseMatrix(const char *file) {

    std::ifstream in(file);

    in >> H >> W;

    V.push_back(0);
    C.push_back(0);
    N.push_back(SPARSE_END);

    F.assign(H+1, SPARSE_END);

    double x;
    size_t c;

    for (size_t r = 1; r <= H; ++r)
    {
        while (true)
        {
            in >> c;

            if (c != SPARSE_END) {
                in >> x;
                add(r, c, x);
            } else {
                break;
            }
        }
    }

    in.close();
}


void SparseMatrix::set(size_t row, size_t col, double value)
{

    if ((row > H) || (col > W)) return;

    size_t j, q, p;

    for (p = SPARSE_END, q = F[row]; q != SPARSE_END; p = q, q = N[q])
    {
        j = C[q];
        if (j < col) continue;
        if (j == col) //set
        {
            V[q] = value;
            return;
        }
        if (j > col) //add
        {
            break;
        }
    }
    V.push_back(value);
    N.push_back(q);
    C.push_back(col);
    if (p == SPARSE_END)
        F[row] = V.size() - 1;
    else
        N[p] = V.size() - 1;
}


double SparseMatrix::get(size_t row, size_t col) const
{
    if ((row > H) || (col > W)) return 0;

    size_t j;

    for (size_t q = F[row]; q != SPARSE_END; q = N[q])
    {
        j = C[q];
        if (j < col) continue;
        if (j == col) return V[q];
        if (j > col) return 0;
    }
    return 0;
}


void SparseMatrix::add(size_t row, size_t col, double value)
{
    size_t q = N.size();
    C.push_back(col);
    V.push_back(value);
    N.push_back(SPARSE_END);
    if (F[row] == SPARSE_END)
        F[row] = q;
    else
        N[q-1] = q;
}


void SparseMatrix::swapRow(size_t r1, size_t r2)
{
    size_t t = F[r1]; F[r1] = F[r2]; F[r2] = t;
}


void SparseMatrix::swapCol(size_t c1, size_t c2)
{

    for (size_t i = 1; i <= H; ++i)
    {
        bool f1 = false, f2 = false;
        size_t q1 = F[i], q2 = F[i];
        size_t p1 = SPARSE_END, p2 = SPARSE_END;
        size_t jmax = std::max(c1,c2);
        for (size_t q = F[i], p = SPARSE_END; q != SPARSE_END; p = q, q = N[q])
        {
            size_t j = C[q];
            if (j == c1) { q1 = q; p1 = p; f1 = true; }
            if (!f1 && (j < c1)) { p1 = q; q1 = N[p1];}

            if (j == c2) { q2 = q; p2 = p; f2 = true; }
            if (!f2 && (j < c2)) { p2 = q; q2 = N[p2];}

            if ((f1 && f2) || (j > jmax)) break;
        }

        if (f1 && f2)
        {
            double t = V[q1];
            V[q1] = V[q2];
            V[q2] = t;
        }
        if (f1 && !f2)
        {
            if ((q1 == q2) || (p2 == q1))
            {
                C[q1] = c2;
            }
            else
            {
                if (p1 == SPARSE_END)
                    F[i] = N[q1];
                else
                    N[p1] = N[q1];

                if (p2 == SPARSE_END)
                    F[i] = q1;
                else
                    N[p2] = q1;
                N[q1] = q2;
                C[q1] = c2;
            }

        }
        if (!f1 && f2)
        {
            if ((q1 == q2) || (p1 == q2))
            {
                C[q2] = c1;
            }
            else
            {
                if (p2 == SPARSE_END)
                    F[i] = N[q2];
                else
                    N[p2] = N[q2];
                if (p1 == SPARSE_END)
                    F[i] = q2;
                else
                    N[p1] = q2;
                N[q2] = q1;
                C[q2] = c1;
            }
        }
    }
}


void SparseMatrix::print() const
{
    for (size_t i = 1; i <= H; ++i)
    {
        size_t j = 1, jj;
        for (size_t q = F[i]; q != SPARSE_END; q = N[q])
        {
            jj = C[q];
            while (j < jj) { std::cout << "0" << '\t'; ++j; }
            std::cout << V[q] << '\t'; ++j;
        }
        while (j <= W) { std::cout << "0" << '\t'; ++j; }
        std::cout << ";\n";
    }
}


void SparseMatrix::save2file(const char *file) const
{
    std::ofstream out(file);

    out << H << " " << W << '\n';

    for (size_t i = 1; i <= H; ++i)
    {
        for (size_t q = F[i]; q != SPARSE_END; q = N[q])
        {
            out << C[q] << " " << V[q] << '\t';
        }
        out << "0\n";
    }
}


void SparseMatrix::save2fileold(const char *file) const
{
    std::ofstream out(file);

    out << H << '\n';

    for (size_t i = 1; i <= H; ++i)
    {
        for (size_t q = F[i]; q != SPARSE_END; q = N[q])
        {
            out << C[q]-1 << " " << V[q] << '\t';
        }
        out << "-1\n";
    }
}


Matrix operator *(const SparseMatrix &S, const Matrix &M)
{
    if (S.W != M.H) {} // надо бы бросить исключение

    Matrix OM(S.H,M.W);

#pragma omp parallel for
    for (size_t i = 1; i <= OM.H; ++i)
    {
        size_t h = OM.W*(i-1);
        for (size_t q = S.F[i]; q != SPARSE_END; q = S.N[q])
        {
            size_t h1 = M.W*(S.C[q]-1);
            double v = S.V[q];
            for (size_t k = 1; k <= OM.W; ++k)
            {
                //OM[i,k] += S[i,j]*M[j,k]
                OM.M[h+k] += v * M.M[h1+k];
            }
        }
    }

    return OM;
}


Matrix operator *(const Matrix &M, const SparseMatrix &S)
{
    if (M.W != S.H) {} // надо бы бросить исключение

    Matrix OM(M.H,S.W);

    for (size_t i = 1; i <= M.W; ++i)
    {
        for (size_t q = S.F[i]; q != SPARSE_END; q = S.N[q])
        {
            size_t j = S.C[q];
            double v = S.V[q];
#pragma omp parallel for
            for (size_t k = 1; k <= OM.H; ++k)
            {
                //OM[k,j] += M[k,i]*S[i,j]
                OM.M[OM.W*(k-1)+j] += M.M[M.W*(k-1)+i] * v;
            }
        }
    }

    return OM;
}


Vector operator *(const SparseMatrix &S, const Vector &V)
{
    if (S.W != V.H) {} // надо бы бросить исключение

    Vector OV(S.H);

#pragma omp parallel for
    for (size_t i = 1; i <= S.H; ++i)
        for (size_t q = S.F[i]; q != SPARSE_END; q = S.N[q])
            OV.V[i] += S.V[q] * V.V[S.C[q]];

    return OV;
}


Matrix operator +(const SparseMatrix &S, const Matrix &M)
{
    if ((S.H != M.H) || (S.W != M.W)) {} // надо бы бросить исключение

    Matrix OM(M);

#pragma omp parallel for
    for (size_t i = 1; i <= S.H; ++i)
    {
        size_t h = OM.W*(i-1);
        for (size_t q = S.F[i]; q != SPARSE_END; q = S.N[q])
            OM.M[h+S.C[q]] += S.V[q];
    }

    return OM;
}


Matrix operator +(const Matrix &M, const SparseMatrix &S)
{
    return S + M;
}


Matrix& operator +=(Matrix &M, const SparseMatrix &S)
{
    if ((S.H != M.H) || (S.W != M.W)) {} // надо бы бросить исключение

#pragma omp parallel for
    for (size_t i = 1; i <= S.H; ++i)
    {
        size_t h = M.W*(i-1);
        for (size_t q = S.F[i]; q != SPARSE_END; q = S.N[q])
            M.M[h+S.C[q]] += S.V[q];
    }

    return M;
}



Matrix operator -(const SparseMatrix &S, const Matrix &M)
{
    if ((S.H != M.H) || (S.W != M.W)) {} // надо бы бросить исключение

    Matrix OM(M);

#pragma omp parallel for
    for (size_t i = 1; i < OM.M.size(); ++i) OM.M[i] = -OM.M[i];

#pragma omp parallel for
    for (size_t i = 1; i <= S.H; ++i)
    {
        size_t h = OM.W*(i-1);
        for (size_t q = S.F[i]; q != SPARSE_END; q = S.N[q])
            OM.M[h+S.C[q]] += S.V[q];
    }

    return OM;
}


Matrix operator -(const Matrix &M, const SparseMatrix &S)
{
    if ((S.H != M.H) || (S.W != M.W)) {} // надо бы бросить исключение

    Matrix OM(M);

#pragma omp parallel for
    for (size_t i = 1; i <= S.H; ++i)
    {
        size_t h = OM.W*(i-1);
        for (size_t q = S.F[i]; q != SPARSE_END; q = S.N[q])
            OM.M[h+S.C[q]] -= S.V[q];
    }

    return OM;
}


#endif // SPARSEMATRIX_H
