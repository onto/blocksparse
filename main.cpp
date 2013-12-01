#include <iostream>

#include "sparsematrix.h"
#include "matrixoperations.h"
#include "decompositor.h"
#include <cstdlib>
#include <time.h>

using namespace std;

int main(int argc, char *argv[3])
{
    if (argc < 3)
    {
        cout << "Need type and threads number" << endl;
        return 0;
    }

    int type = atoi(argv[1]);
    int threads = atoi(argv[2]);

    if (type == 1)
    {
        SparseMatrix M("matrix.txt");

        Vector B("vector.txt");

        ofstream out("st.txt",ios_base::app);

#ifdef _OPENMP
        omp_set_num_threads(threads);
        double t = omp_get_wtime();
        double t1 = t;
#else
        time_t t = clock();
        time_t t1 = t;
#endif

        LUPS T;
        MatrixOperations::LUTriang(M, T);

#ifdef _OPENMP
        cout << "Время на разложение " << omp_get_wtime() - t << endl;
        t = omp_get_wtime();
#else
        cout << "Время на разложение " << double((clock()-t))/CLOCKS_PER_SEC << endl;
        t = clock();
#endif

        Vector X;
        MatrixOperations::Solve(T, B, X);

#ifdef _OPENMP
        cout << "Время на решение " << omp_get_wtime() - t << endl;
        cout << "Общее время " << omp_get_wtime() - t1 << endl;
        //out << omp_get_wtime() - t1 << endl;
#else
        cout << "Время на решение " << double((clock()-t))/CLOCKS_PER_SEC << endl;
        cout << "Общее время " << double((clock()-t1))/CLOCKS_PER_SEC << endl;
        //out << double((clock()-t1))/CLOCKS_PER_SEC << endl;
#endif

//        X.print();
        cout << X.V[1] << endl;
        cout << X.V[2] << endl;
    }
    else if (type == 2)
    {
        BlockSparseMatrix M1("matrix_b.txt", BlockSparseMatrix::BlockMatrixInputType);
        Vector B("vector.txt");

        size_t S = 0;
        for (size_t i = 0; i < M1.A.size(); ++i)
        {
            S += M1.A[i].V.size()-1;
            S += M1.B[i].V.size()-1;
            S += M1.C[i].V.size()-1;
        }
        S += M1.Q.V.size()-1;

        //cout << S << endl;

#ifdef _OPENMP
        omp_set_num_threads(threads);
        double t = omp_get_wtime();
        double t1 = t;
        ofstream out;
        if (threads > 1)
            out.open("ot.txt",ios_base::app);
        else
            out.open("bt.txt",ios_base::app);
#else
        time_t t = clock();
        time_t t1 = t;
        ofstream out("bt.txt",ios_base::app);
#endif

        FactorizedBlockSparseMatrix FM;
        MatrixOperations::BlockMatrixFactorization(M1, FM);

#ifdef _OPENMP
        //cout << "Время на разложение " << omp_get_wtime() - t << endl;
        t = omp_get_wtime();
#else
        //cout << "Время на разложение " << double((clock()-t))/CLOCKS_PER_SEC << endl;
        t = clock();
#endif

        Vector X;
        MatrixOperations::Solve(FM, B, X);

#ifdef _OPENMP
        //cout << "Время на решение " << omp_get_wtime() - t << endl;
       // cout << "Общее время " << omp_get_wtime() - t1 << endl;
        out << omp_get_wtime() - t1 << endl;
#else
        //cout << "Время на решение " << double((clock()-t))/CLOCKS_PER_SEC << endl;
        //cout << "Общее время " << double((clock()-t1))/CLOCKS_PER_SEC << endl;
        out << double((clock()-t1))/CLOCKS_PER_SEC << endl;
#endif

        cout << X.V[1] << endl;
        cout << X.V[2] << endl;
    }
    else if (type == 3)
    {

#ifdef _OPENMP
        omp_set_num_threads(threads);
        double t = omp_get_wtime();
        double t1 = t;
#else
        time_t t = clock();
        time_t t1 = t;
#endif

        BlockSparseMatrix M("matrix.txt", BlockSparseMatrix::SparseMatrixInputType);

#ifdef _OPENMP
        cout << "Время на декомпозицию " << omp_get_wtime() - t << endl;
        t = omp_get_wtime();
#else
        cout << "Время на декомпозицию " << double((clock()-t))/CLOCKS_PER_SEC << endl;
        t = clock();
#endif

        FactorizedBlockSparseMatrix FM;
        MatrixOperations::BlockMatrixFactorization(M, FM);

#ifdef _OPENMP
        cout << "Время на разложение " << omp_get_wtime() - t << endl;
        t = omp_get_wtime();
#else
        cout << "Время на разложение " << double((clock()-t))/CLOCKS_PER_SEC << endl;
        t = clock();
#endif

        Vector B("vector.txt");
        B.permute(M.D.Pr);

        Vector X;

        MatrixOperations::Solve(FM, B, X);

        X.permute(M.D.Pct);

#ifdef _OPENMP
        cout << "Время на решение " << omp_get_wtime() - t << endl;
        cout << "Общее время " << omp_get_wtime() - t1 << endl;
        //out << omp_get_wtime() - t1 << endl;
#else
        cout << "Время на решение " << double((clock()-t))/CLOCKS_PER_SEC << endl;
        cout << "Общее время " << double((clock()-t1))/CLOCKS_PER_SEC << endl;
        //out << double((clock()-t1))/CLOCKS_PER_SEC << endl;
#endif

        //X.print();
        cout << X.V[1] << endl;
        cout << X.V[2] << endl;

    }
    else if (type == 4)
    {

        Matrix M("matrix.txt");

#ifdef _OPENMP
        omp_set_num_threads(threads);
        double t = omp_get_wtime();
        double t1 = t;
#else
        time_t t = clock();
        time_t t1 = t;
#endif

        LUPM T;
        MatrixOperations::LUTriang(M, T);

#ifdef _OPENMP
        cout << "Время на разложение " << omp_get_wtime() - t << endl;
        t = omp_get_wtime();
#else
        cout << "Время на разложение " << double((clock()-t))/CLOCKS_PER_SEC << endl;
        t = clock();
#endif

        Vector B("vector.txt");
        Vector X;
        MatrixOperations::Solve(T, B, X);

#ifdef _OPENMP
        cout << "Время на решение " << omp_get_wtime() - t << endl;
        cout << "Общее время " << omp_get_wtime() - t1 << endl;
#else
        cout << "Время на решение " << double((clock()-t))/CLOCKS_PER_SEC << endl;
        cout << "Общее время " << double((clock()-t1))/CLOCKS_PER_SEC << endl;
#endif

//        X.print();
        cout << X.V[1] << endl;
        cout << X.V[2] << endl;
    }
    else if (type == 5)
    {

#ifdef _OPENMP
        omp_set_num_threads(threads);
        double t = omp_get_wtime();
        double t1 = t;
#else
        time_t t = clock();
        time_t t1 = t;
#endif

        BlockSparseMatrix M("matrix.txt", BlockSparseMatrix::SparseMatrixUnsimmetricInputType);

#ifdef _OPENMP
        cout << "Время на декомпозицию " << omp_get_wtime() - t << endl;
        t = omp_get_wtime();
#else
        cout << "Время на декомпозицию " << double((clock()-t))/CLOCKS_PER_SEC << endl;
        t = clock();
#endif

        FactorizedBlockSparseMatrix FM;
        MatrixOperations::BlockMatrixFactorization(M, FM);

#ifdef _OPENMP
        cout << "Время на разложение " << omp_get_wtime() - t << endl;
        t = omp_get_wtime();
#else
        cout << "Время на разложение " << double((clock()-t))/CLOCKS_PER_SEC << endl;
        t = clock();
#endif

        Vector B("vector.txt");
        B.permute(M.D.Pr);

        Vector X;

        MatrixOperations::Solve(FM, B, X);

        X.permute(M.D.Pct);

#ifdef _OPENMP
        cout << "Время на решение " << omp_get_wtime() - t << endl;
        cout << "Общее время " << omp_get_wtime() - t1 << endl;
        //out << omp_get_wtime() - t1 << endl;
#else
        cout << "Время на решение " << double((clock()-t))/CLOCKS_PER_SEC << endl;
        cout << "Общее время " << double((clock()-t1))/CLOCKS_PER_SEC << endl;
        //out << double((clock()-t1))/CLOCKS_PER_SEC << endl;
#endif

        //X.print();
        cout << X.V[1] << endl;
        cout << X.V[2] << endl;

    }

//    SparseMatrix M("A0.txt");
//    Matrix Minv;
//    MatrixOperations::Inverse(M, Minv);
//    Minv.print();

    return 0;
}
