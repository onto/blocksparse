#ifndef LUCONTAINER_H
#define LUCONTAINER_H

#include <vector>

#include "matrix.h"
#include "sparsematrix.h"
#include "decompositor.h"

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

    void save2fileL(const char *file)
    {
        std::ofstream out(file);

        out << U.H << " " << U.W << '\n';

        for (size_t i = 1; i <= U.H; ++i)
        {
            for (size_t q = LF[i]; q != SPARSE_END; q = U.N[q])
            {
                out << U.C[q] << " " << U.V[q] << '\t';
            }
            out << "0\n";
        }
        out.close();
    }

    void save2fileU(const char *file)
    {
        U.save2file(file);
    }
};

struct LUPM
{
    Matrix M;
    std::vector<size_t> P;
    std::vector<size_t> Pt;

//    void printL() const
//    {
//        for (size_t i = 1; i <= M.H; ++i)
//        {
//            for (size_t j = 1; j < i; ++j)
//                std::cout << M(i,j) << '\t';
//            for (size_t j = i; j <= M.H; ++j)
//                std::cout << "*\t";
//            std::cout << '\n';
//        }
//    }

//    void printU() const
//    {
//        for (size_t i = 1; i <= M.H; ++i)
//        {
//            for (size_t j = 1; j < i; ++j)
//                std::cout << "*\t";
//            for (size_t j = i; j <= M.H; ++j)
//                std::cout << M(i,j) << '\t';
//            std::cout << '\n';
//        }
//    }
};

struct FactorizedBBDSparseMatrix
{
    std::vector< LUPS > Alu;
    std::vector< Matrix > Bh;
    std::vector< SparseMatrix > C;
    LUPM H;
    BBDStruct BBDS;

    size_t N, Nb;
};

#endif // LUCONTAINER_H
