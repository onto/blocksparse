#ifndef BLOCKMATRIX_H
#define BLOCKMATRIX_H

#include <vector>
#include <iostream>
#include <fstream>

#include "sparsematrix.h"

struct BlockSparseMatrix
{
    std::vector< SparseMatrix > A;   // A_i - r_i x r_i
    std::vector< SparseMatrix > B;   // B_i - r_i x r_q
    std::vector< SparseMatrix > C;   // C_i - r_q x r_i
    SparseMatrix Q;                  // Q - r_q x r_q
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

        double x;
        size_t j;

        size_t q = 0, q1 = N - R.back();
        for (size_t k = 0; k < R.size()-1; ++k)
        {
            A.push_back(SparseMatrix(R[k], R[k]));
            B.push_back(SparseMatrix(R[k], R.back()));
            C.push_back(SparseMatrix(R.back(), R[k]));
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

        Q = SparseMatrix(R.back(), R.back());

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

#endif // BLOCKMATRIX_H
