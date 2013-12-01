#ifndef LUCONTAINER_H
#define LUCONTAINER_H

#include <vector>

#include "matrix.h"
#include "sparsematrix.h"

struct LUPS
{
    SparseMatrix U;
    std::vector<size_t> LF;
    std::vector<size_t> P;
    std::vector<size_t> Pt;

    void printL() const
    {
        for (size_t i = 1; i <= U.H; ++i)
        {
            size_t j = 1, jj;
            for (size_t q = LF[i]; q != SPARSE_END; q = U.N[q])
            {
                jj = U.C[q];
                while (j < jj) { std::cout << "0" << '\t'; ++j; }
                std::cout << U.V[q] << '\t'; ++j;
            }
            while (j <= U.W) { std::cout << "0" << '\t'; ++j; }
            std::cout << ";\n";
        }
    }

    void printU() const
    {
        U.print();
    }
};

struct LUPM
{
    Matrix M;
    std::vector<size_t> P;
    std::vector<size_t> Pt;

    void printL() const
    {
        for (size_t i = 1; i <= M.H; ++i)
        {
            for (size_t j = 1; j < i; ++j)
                std::cout << M.get(i,j) << '\t';
            for (size_t j = i; j <= M.H; ++j)
                std::cout << "*\t";
            std::cout << '\n';
        }
    }

    void printU() const
    {
        for (size_t i = 1; i <= M.H; ++i)
        {
            for (size_t j = 1; j < i; ++j)
                std::cout << "*\t";
            for (size_t j = i; j <= M.H; ++j)
                std::cout << M.get(i,j) << '\t';
            std::cout << '\n';
        }
    }
};

struct FactorizedBlockSparseMatrix
{
    std::vector< LUPS > Alu;
    std::vector< Matrix > Bh;
    std::vector< SparseMatrix > C;
    LUPM H;
    std::vector<size_t> R;
};

#endif // LUCONTAINER_H
