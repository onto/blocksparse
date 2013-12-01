#ifndef BLOCKMATRIX_H
#define BLOCKMATRIX_H

#include <vector>
#include <iostream>
#include <fstream>

#include "sparsematrix.h"
#include "decompositor.h"
#include "udecompositor.h"

struct BlockSparseMatrix
{
    std::vector< SparseMatrix > A;   // A_i - r_i x r_i
    std::vector< SparseMatrix > B;   // B_i - r_i x r_q
    std::vector< SparseMatrix > C;   // C_i - r_q x r_i
    SparseMatrix Q;                  // Q - r_q x r_q
    std::vector<size_t> R;           // r_1--r_q

    SimpleMatrixDecompositor D;

    bool unsimmetric;

    enum InputType
    {
        BlockMatrixInputType = 0,
        SparseMatrixInputType,
        SparseMatrixUnsimmetricInputType
    };

    BlockSparseMatrix() {}

    BlockSparseMatrix(const char *file, InputType type)
    {
        switch (type) {
        case BlockMatrixInputType:
            readBlockMatrix(file);
            break;
        case SparseMatrixInputType:
            unsimmetric = false;
            readSparseMatrixAndDecompose(file);
            break;
        case SparseMatrixUnsimmetricInputType:
            unsimmetric = true;
            readUnsimmetricSparseMatrixAndDecompose(file);
            break;
        default:
            break;
        }
    }

private:
    void readBlockMatrix(const char *file)
    {
        std::ifstream in(file);

        R.clear();

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

        size_t q = 0, qn = 0, q1 = N - R.back();
        for (size_t k = 0; k < R.size()-1; ++k)
        {
            q = qn;
            qn += R[k];

            A.push_back(SparseMatrix(R[k], R[k]));
            B.push_back(SparseMatrix(R[k], R.back()));
            C.push_back(SparseMatrix(R.back(), R[k]));
            for (size_t i = 1; i <= R[k]; ++i)
            {
                in >> j;
                while (j != SPARSE_END)
                {
                    in >> x;
                    if (j <= qn)
                        A[k].add(i,j-q,x);
                    else
                        B[k].add(i,j-q1,x);
                    in >> j;
                }
            }
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

    void readSparseMatrixAndDecompose(const char *file)
    {
        MatrixDecompositor Dt(file);
        D = Dt;
        R = Dt.R;

        //Теперь считываем, производя перестановку
        SparseMatrix S(file);
        S.permute(D.Pr);

        S.save2file("matrix.perm.txt");

//        std::ofstream out("matrix.temp.txt");

//        for (size_t i = 0; i < R.size(); ++i)
//            out << R[i] << " ";
//        out <<"0\n";

//        for (size_t i = 1; i <= S.H; ++i)
//        {
//            for (size_t q = S.F[i]; q != SPARSE_END; q = S.N[q])
//            {
//                out << S.C[q] << " " << S.V[q] << '\t';
//            }
//            out << "0\n";
//        }
//        out.close();

//        readBlockMatrix("matrix.temp.txt");

        //Теперь надо заполнить A, B, C, Q
        size_t qn = 0, qp = 0, c, qb = S.H - R.back();
        for (size_t k = 0; k < R.size()-1; ++k)
        {
            qp = qn;
            qn += R[k];

            A.push_back(SparseMatrix(R[k], R[k]));
            B.push_back(SparseMatrix(R[k], R.back()));
            C.push_back(SparseMatrix(R.back(), R[k]));
            for (size_t i = 1; i <= R[k]; ++i)
            {
                for (size_t q = S.F[qp+i]; q != SPARSE_END; q = S.N[q])
                {
                    c = S.C[q];
                    if (c <= qn)
                        A[k].add(i, c-qp, S.V[q]);
                    else
                        B[k].add(i, c-qb, S.V[q]);
                }
            }
        }

        Q = SparseMatrix(R.back(),R.back());

        for (size_t i = 1; i <= R.back(); ++i)
        {
            qp = 0; qn = R[0];
            size_t k = 0, q;
            for (q = S.F[qb+i]; q != SPARSE_END; q = S.N[q])
            {
                c = S.C[q];
                while (c > qn)
                {
                    qp = qn;
                    ++k;
                    qn += R[k];
                }
                if (k == R.size()-1) break;

                C[k].add(i, c-qp, S.V[q]);
            }

            for (; q != SPARSE_END; q = S.N[q])
            {
                Q.add(i, S.C[q] - qb, S.V[q]);
            }
        }
    }

    void readUnsimmetricSparseMatrixAndDecompose(const char *file)
    {
        MatrixUDecompositor Dt(file);
        D = Dt;
        R = Dt.R;

        //Теперь считываем, производя перестановку
        SparseMatrix S(file);
        S.permute(D.Pc, D.Pc);

        S.save2file("matrix.perm.txt");

//        std::ofstream out("matrix.temp.txt");

//        for (size_t i = 0; i < R.size(); ++i)
//            out << R[i] << " ";
//        out <<"0\n";

//        for (size_t i = 1; i <= S.H; ++i)
//        {
//            for (size_t q = S.F[i]; q != SPARSE_END; q = S.N[q])
//            {
//                out << S.C[q] << " " << S.V[q] << '\t';
//            }
//            out << "0\n";
//        }
//        out.close();

//        readBlockMatrix("matrix.temp.txt");

        //Теперь надо заполнить A, B, C, Q
        size_t qn = 0, qp = 0, c, qb = S.H - R.back();
        for (size_t k = 0; k < R.size()-1; ++k)
        {
            qp = qn;
            qn += R[k];

            A.push_back(SparseMatrix(R[k], R[k]));
            B.push_back(SparseMatrix(R[k], R.back()));
            C.push_back(SparseMatrix(R.back(), R[k]));
            for (size_t i = 1; i <= R[k]; ++i)
            {
                for (size_t q = S.F[qp+i]; q != SPARSE_END; q = S.N[q])
                {
                    c = S.C[q];
                    if (c <= qn)
                        A[k].add(i, c-qp, S.V[q]);
                    else
                        B[k].add(i, c-qb, S.V[q]);
                }
            }
        }

        Q = SparseMatrix(R.back(),R.back());

        for (size_t i = 1; i <= R.back(); ++i)
        {
            qp = 0; qn = R[0];
            size_t k = 0, q;
            for (q = S.F[qb+i]; q != SPARSE_END; q = S.N[q])
            {
                c = S.C[q];
                while (c > qn)
                {
                    qp = qn;
                    ++k;
                    qn += R[k];
                }
                if (k == R.size()-1) break;

                C[k].add(i, c-qp, S.V[q]);
            }

            for (; q != SPARSE_END; q = S.N[q])
            {
                Q.add(i, S.C[q] - qb, S.V[q]);
            }
        }
    }
};

#endif // BLOCKMATRIX_H
