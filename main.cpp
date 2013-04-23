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

        LUPS T;

        time_t t = clock();

        MatrixOperations::LUTriang(M, T);

        cout << "Время на разложение " << double((clock()-t))/CLOCKS_PER_SEC << endl;

        Vector X;

        t = clock();
        MatrixOperations::Solve(T, B, X);

        cout << "Время на решение " << double((clock()-t))/CLOCKS_PER_SEC << endl;

        //X.print();
        cout << X.V[1] << endl;
        cout << X.V[2] << endl;
    }
    else
    {
        BlockSparseMatrix M1("matrix_b.txt");
        Vector B("vector.txt");

#ifdef _OPENMP
        omp_set_num_threads(threads);
        double t = omp_get_wtime();
#else
    time_t t = clock();
#endif

        FactorizedBlockSparseMatrix FM;
        MatrixOperations::BlockMatrixFactorization(M1, FM);

#ifdef _OPENMP
        t = omp_get_wtime() - t;
        cout << "Время на разложение " << t << endl;
        t = omp_get_wtime();
#else
        cout << "Время на разложение " << double((clock()-t))/CLOCKS_PER_SEC << endl;
        t = clock();
#endif

        Vector X;
        MatrixOperations::Solve(FM, B, X);

#ifdef _OPENMP
        t = omp_get_wtime() - t;
        cout << "Время на решение " << t << endl;
#else
        cout << "Время на решение " << double((clock()-t))/CLOCKS_PER_SEC << endl;
#endif
        cout << X.V[1] << endl;
        cout << X.V[2] << endl;
    }

//    SparseMatrix M("m6.txt");
//    M.print();

//    Matrix Minv;

//    MatrixOperations::Inverse(M, Minv);

//    Minv.print();

    return 0;
}
