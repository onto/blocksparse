#ifndef BLOCKMATRIX_H
#define BLOCKMATRIX_H

#include <vector>
#include <iostream>
#include <fstream>

#include <matrixoperations.h>
#include <sparsematrix.h>
#include <matrix.h>

template <typename T> struct LUPM;

template <typename T> struct BlockSparseMatrix
{
    std::vector< SparseMatrix<T> > A;   // A_i - r_i x r_i
    std::vector< SparseMatrix<T> > B;   // B_i - r_i x r_q
    std::vector< SparseMatrix<T> > C;   // C_i - r_q x r_i
    SparseMatrix<T> Q;                  // Q - r_q x r_q
    std::vector<size_t> R;              // r_1--r_q

    BlockSparseMatrix(const char *file)
    {

        std::ifstream in(file);

        size_t r, N = 0;
        while (true)
        {
            in >> r;
            N += r;
            if (r != SPARSE_END)
                R.push_back(r);
            else
                break;
        }

        T x;
        size_t j;

        size_t q = 0, q1 = N - R.back();
        for (size_t k = 0; k < R.size()-1; ++k)
        {
            A.push_back(SparseMatrix<T>(R[k], R[k]));
            B.push_back(SparseMatrix<T>(R[k], R.back()));
            C.push_back(SparseMatrix<T>(R.back(), R[k]));
            for (size_t i = 1; i <= R[k]; ++i)
            {
                in >> j;
                while (j != SPARSE_END)
                {
                    in >> x;
                    if (j <= q+R[k])
                        A[k].add(i,j-q,x);
                    else
                        B[k].add(i,j-q1,x);
                    in >> j;
                }
            }
            q += R[k];
        }

        Q = SparseMatrix<T>(R.back(), R.back());

        for (size_t i = 1; i <= R.back(); ++i)
        {
            q = 0;
            size_t k = 0;
            in >> j;
            while (j != SPARSE_END)
            {
                while (j > q+R[k]) {
                    q += R[k];
                    ++k;
                }
                if (k == R.size()-1) break;
                in >> x;
                C[k].add(i,j-q,x);
                in >> j;
            }

            //in >> j;
            while (j != SPARSE_END)
            {
                in >> x;
                Q.add(i,j-q1,x);
                in >> j;
            }

        }

        in.close();
    }
};

template <typename T> struct FactorizedBlockSparseMatrix
{
    std::vector< Matrix<T> > Ainv;      // Ainv_i = inv(A_i)
    std::vector< SparseMatrix<T> > B;   // B_i - r_i x r_q
    std::vector< SparseMatrix<T> > C;   // C_i - r_q x r_i
    SparseMatrix<T> Q;                  // Q - r_q x r_q
    LUPM<T> H;                          // H = LUTriang(Q - sum(C_i * Ainv_i * B_i))
    std::vector<size_t> R;              // r_1--r_q
};

#endif // BLOCKMATRIX_H
