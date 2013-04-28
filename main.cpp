#include <iostream>

#include "sparsematrix.h"
#include "matrixoperations.h"
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
        SparseMatrix M("matrix_s.txt");

        Vector B("vector.txt");

        //M.print(); B.print();

        cout << M.N.size()-1 << endl;

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
#else
        cout << "Время на решение " << double((clock()-t))/CLOCKS_PER_SEC << endl;
        cout << "Общее время " << double((clock()-t1))/CLOCKS_PER_SEC << endl;
#endif

        //X.print();
        cout << X.V[1] << endl;
        cout << X.V[2] << endl;
    }
    else
    {
        BlockSparseMatrix M1("matrix_b.txt");
        Vector B("vector.txt");

        size_t S = 0;
        for (size_t i = 0; i < M1.A.size(); ++i)
        {
            S += M1.A[i].V.size()-1;
            S += M1.B[i].V.size()-1;
            S += M1.C[i].V.size()-1;
        }
        S += M1.Q.V.size()-1;

        cout << S << endl;

#ifdef _OPENMP
        omp_set_num_threads(threads);
        double t = omp_get_wtime();
        double t1 = t;
#else
        time_t t = clock();
        time_t t1 = t;
#endif

        FactorizedBlockSparseMatrix FM;
        MatrixOperations::BlockMatrixFactorization(M1, FM);

#ifdef _OPENMP
        cout << "Время на разложение " << omp_get_wtime() - t << endl;
        t = omp_get_wtime();
#else
        cout << "Время на разложение " << double((clock()-t))/CLOCKS_PER_SEC << endl;
        t = clock();
#endif

        Vector X;
        MatrixOperations::Solve(FM, B, X);

#ifdef _OPENMP
        cout << "Время на решение " << omp_get_wtime() - t << endl;
        cout << "Общее время " << omp_get_wtime() - t1 << endl;
#else
        cout << "Время на решение " << double((clock()-t))/CLOCKS_PER_SEC << endl;
        cout << "Общее время " << double((clock()-t1))/CLOCKS_PER_SEC << endl;
#endif

        cout << X.V[1] << endl;
        cout << X.V[2] << endl;
    }

//    SparseMatrix M("matrix_s.txt");
//    //M.print();

//    Matrix Minv;

//    time_t t = clock();
//    MatrixOperations::Inverse(M, Minv);
//    cout << "Время на решение " << double((clock()-t))/CLOCKS_PER_SEC << endl;

    //Minv.print();

    return 0;
}
