#ifndef BLOCKMATRIX_H
#define BLOCKMATRIX_H

#include <vector>
#include <iostream>
#include <fstream>

#include "sparsematrix.h"
#include "decompositor.h"

struct BBDSparseMatrix
{
    std::vector< SparseMatrix > A;   // A_i - r_i x r_i
    std::vector< SparseMatrix > B;   // B_i - r_i x r_q
    std::vector< SparseMatrix > C;   // C_i - r_q x r_i
    SparseMatrix D;                  // D   - r_q x r_q
    BBDStruct BBDS;                  // Структура, хранящая информацию о
                                     // перестановке исходной матрицы и др.

    size_t N, Nb; //Общий размер, количество блоков

    struct ComparePairByFirst
    {
        bool operator()(const std::pair<size_t, size_t> &p1, const std::pair<size_t, size_t> &p2)
        {
            return p1.first < p2.first;
        }
    } comparator;

    BBDSparseMatrix(SparseMatrix M)
    {
        N = M.H;

        //Определяем перестановку и структуру
//        MatrixBBDPreordering::BBDDecompose(M, BBDS, SV, 0.01875*N, 0.025*N, 0.025*N);
//        showBlocksInfo();
//        MatrixBBDPreordering::BBDDecompose(M, BBDS, SV, 0.0375*N, 0.05*N, 0.05*N);
//        showBlocksInfo();
        MatrixBBDPreordering::BBDDecompose(M, BBDS, SV, 0.075*N, 0.1*N, 0.5*N);
        printBlocksInfo();

//        MatrixBBDPreordering::BBDDecompose(M, BBDS, SVI, 0.01875*N, 0.025*N, 0.025*N);
//        showBlocksInfo();
//        MatrixBBDPreordering::BBDDecompose(M, BBDS, SVI, 0.0375*N, 0.05*N, 0.05*N);
//        showBlocksInfo();
        MatrixBBDPreordering::BBDDecompose(M, BBDS, SVI, 0.075*N, 0.1*N, 0.5*N);
        printBlocksInfo();

        //Переставляем столбцы исходной матрицы
        for (size_t i = 1; i <= N; ++i)
        {
            std::vector< std::pair<size_t, size_t> > CQ;

            for (size_t q = M.F[i]; q != SPARSE_END; q = M.N[q])
            {
                CQ.push_back( std::pair<size_t, size_t>(BBDS.Pct[M.C[q]], q) );
            }

            std::sort(CQ.begin(), CQ.end(), comparator);

            //Меняем порядок
            M.F[i] = CQ[0].second; //Новое начало строки

            size_t cqend = CQ.size()-1;
            for (size_t k = 0; k < cqend; ++k)
            {
                M.C[CQ[k].second] = CQ[k].first;
                M.N[CQ[k].second] = CQ[k+1].second;
            }
            //Последнему элементу строки в A.N назначается SPARSE_END
            M.C[CQ[cqend].second] = CQ[cqend].first;
            M.N[CQ[cqend].second] = SPARSE_END;
        }

        //Переставляем строки
        VectSizet Ft(N+1);
        for (size_t i = 1; i <= N; ++i) Ft[i] = M.F[BBDS.Pr[i]];
        M.F.swap(Ft);

        //M.save2file("matrix.perm.txt");
        //saveBlocks2File("matrix.blocks.txt");
        //saveInterfaceSizes2File("matrix.is.txt");

        //Осталось заполнить блоки

        Nb = BBDS.R.size() - 2; //Количество диагональных блоков
        size_t adim = BBDS.R[Nb]; //Суммарная размерность блоков A_i
        size_t ddim = N - adim; //Размерность блока D


        //Заполняем A_i, B_i
        for (size_t ib = 0; ib < Nb; ++ib)
        {
            size_t b0 = BBDS.R[ib], bb = b0+1, be = BBDS.R[ib+1], j;
            size_t bdim = be - b0;

            A.push_back( SparseMatrix(bdim, bdim) );
            B.push_back( SparseMatrix(bdim, ddim) );
            C.push_back( SparseMatrix(ddim, bdim) );

            for (size_t i = bb; i <= be; ++i)
            {
                for (size_t q = M.F[i]; q != SPARSE_END; q = M.N[q])
                {
                    if ((j = M.C[q]) <= be)
                    {
                        A[ib].add(i-b0, j-b0, M.V[q]);
                    }
                    else
                    {
                        B[ib].add(i-b0, j-adim, M.V[q]);
                    }
                }
            }
        }

        D = SparseMatrix(ddim, ddim);

        //Заполняем C_i, D
        for (size_t i = adim+1; i <= N; ++i)
        {
            size_t ib = 0, j;
            for (size_t q = M.F[i]; q != SPARSE_END; q = M.N[q])
            {
                j = M.C[q];
                while ((ib < Nb) && !((j > BBDS.R[ib]) && (j <= BBDS.R[ib+1]))) { ++ib; }

                size_t b0 = BBDS.R[ib];

                if (ib == Nb)
                {
                    D.add(i-adim, j-adim, M.V[q]);
                }
                else
                {
                    C[ib].add(i-adim, j-b0, M.V[q]);
                }
            }
        }
    }

    void printBlocks()
    {
        for (size_t i = 1; i < BBDS.R.size(); ++i)
            std::cout << BBDS.R[i]-BBDS.R[i-1] << " ";
        std::cout << std::endl;
    }

    void saveBlocks2File(const char *file)
    {
        std::ofstream out(file);

        for (size_t i = 1; i < BBDS.R.size(); ++i)
            out << BBDS.R[i]-BBDS.R[i-1] << " ";
        out << std::endl;

        out.close();
    }

    void saveInterfaceSizes2File(const char *file)
    {
        std::ofstream out(file);

        for (size_t i = 0; i < BBDS.Is.size(); ++i)
        {
            out << BBDS.Is[i] << " ";
        }
        out << std::endl;

        out.close();
    }

    void printBlocksInfo()
    {
        size_t maxBlock = 0;

        for (size_t i = 1; i < BBDS.R.size()-1; ++i)
            maxBlock = std::max(BBDS.R[i]-BBDS.R[i-1], maxBlock);

        std::cout << "Максимальный блок A: " << maxBlock << std::endl;
        std::cout << "Размер блока D: " << BBDS.R[BBDS.R.size()-1]-BBDS.R[BBDS.R.size()-2] << std::endl;
        std::cout << std::endl;
    }
};

#endif // BLOCKMATRIX_H
