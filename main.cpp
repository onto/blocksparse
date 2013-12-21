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

    ofstream out("res.csv",ios_base::app);

    if (type == 1)
    {
        SparseMatrix M("matrix.txt");

#ifdef _OPENMP
        omp_set_num_threads(threads);
        double t = omp_get_wtime();
        double t1 = t;
        double tdec, tsol, tres;
#else
        time_t t = clock();
        time_t t1 = t;
        time_t tdec, tsol, tres;
#endif

        LUPS T;
        MatrixOperations::LUTriang(M, T);

#ifdef _OPENMP
        tdec = omp_get_wtime() - t;
        cout << "Время на разложение " << tdec << endl;
        t = omp_get_wtime();
#else
        tdec = double((clock()-t))/CLOCKS_PER_SEC
        cout << "Время на разложение " << tdec << endl;
        t = clock();
#endif

        Vector B("vector.txt");
        Vector X;
        MatrixOperations::Solve(T, B, X);

#ifdef _OPENMP
        tsol = omp_get_wtime() - t;
        tres = omp_get_wtime() - t1;
        cout << "Время на решение " << tsol << endl;
        cout << "Общее время " << tres << endl;
#else
        tsol = double((clock()-t))/CLOCKS_PER_SEC;
        tres = double((clock()-t1))/CLOCKS_PER_SEC;
        cout << "Время на решение " << tsol << endl;
        cout << "Общее время " << tres << endl;
#endif

        //out << tdec << "; " << tsol << "; " << tres << "; ";

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
        double tper, tdec, tsol, tres;
#else
        time_t t = clock();
        time_t t1 = t;
        time_t tper, tdec, tsol, tres;
#endif

        BlockSparseMatrix M("matrix.txt", BlockSparseMatrix::SparseMatrixInputType);

#ifdef _OPENMP
        tper = omp_get_wtime() - t;
        cout << "Время на декомпозицию " << tper << endl;
        t = omp_get_wtime();
#else
        tper = double((clock()-t))/CLOCKS_PER_SEC;
        cout << "Время на декомпозицию " << tper << endl;
        t = clock();
#endif

        FactorizedBlockSparseMatrix FM;
        MatrixOperations::BlockMatrixFactorization(M, FM);

#ifdef _OPENMP
        tdec = omp_get_wtime() - t;
        cout << "Время на разложение " << tdec << endl;
        t = omp_get_wtime();
#else
        tdec = double((clock()-t))/CLOCKS_PER_SEC
        cout << "Время на разложение " << tdec << endl;
        t = clock();
#endif

        Vector B("vector.txt");
        B.permute(M.D.Pr);

        Vector X;

        MatrixOperations::Solve(FM, B, X);

        X.permute(M.D.Pct);

#ifdef _OPENMP
        tsol = omp_get_wtime() - t;
        tres = omp_get_wtime() - t1;
        cout << "Время на решение " << tsol << endl;
        cout << "Общее время " << tres << endl;
#else
        tsol = double((clock()-t))/CLOCKS_PER_SEC;
        tres = double((clock()-t1))/CLOCKS_PER_SEC;
        cout << "Время на решение " << tsol << endl;
        cout << "Общее время " << tres << endl;
#endif

        //out << tper << "; " << tdec << "; " << tsol << "; " << tres << "; ";

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
