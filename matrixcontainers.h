#ifndef LUCONTAINER_H
#define LUCONTAINER_H

#include <vector>

#include <matrix.h>
#include <sparsematrix.h>

struct LUPS
{
    SparseMatrix U;
    std::vector<size_t> LF;
    std::vector<size_t> P;

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
    std::vector< Matrix > Ainv;      // Ainv_i = inv(A_i)
    std::vector< SparseMatrix > B;   // B_i - r_i x r_q
    std::vector< SparseMatrix > C;   // C_i - r_q x r_i
    SparseMatrix Q;                  // Q - r_q x r_q
    LUPM H;                          // H = LUTriang(Q - sum(C_i * Ainv_i * B_i))
    std::vector<size_t> R;              // r_1--r_q
};


#endif // LUCONTAINER_H
