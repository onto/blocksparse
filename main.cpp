#include <iostream>

#ifndef TIME_LOG
#define TIME_LOG

#ifndef CSV_LOG
#define CSV_LOG
#endif

#endif

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


#ifdef CSV_LOG
    ofstream out("res.csv",ios_base::app);
#endif

    if (type == 1) //Решение с разреженной матрицей
    {
        SparseMatrix M("matrix.txt");

#ifdef TIME_LOG
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
#endif

        cout << "Было ненулевых " << M.C.size()-1 << endl;

        LUPS T;
        MatrixOperations::LUTriang(M, T);

        T.save2fileL("spL.txt");
        T.save2fileU("spU.txt");

        cout << "Стало ненулевых " << T.U.C.size()-1 << endl;

#ifdef TIME_LOG
#ifdef _OPENMP
        tdec = omp_get_wtime() - t;
        cout << "Время на разложение " << tdec << endl;
        t = omp_get_wtime();
#else
        tdec = double((clock()-t))/CLOCKS_PER_SEC;
        cout << "Время на разложение " << tdec << endl;
        t = clock();
#endif
#endif

        Vector B("vector.txt"), X;
        MatrixOperations::Solve(T, B, X);

#ifdef TIME_LOG
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

#ifdef CSV_LOG
        out << /*tdec << "; " << tsol << "; " <<*/ tres << "; " << endl;
#endif
#endif

        //X.save2file("solution.txt");

        cout << "Ошибка: " << MatrixOperations::MaxResidual(M, B, X) << endl << endl;

        cout << X.V[1] << endl;
        cout << X.V[2] << endl;
    }
    else if (type == 2) //Решение приведением к БДО форме
    {

        SparseMatrix M("matrix.txt");
        Vector B("vector.txt"), X;

#ifdef TIME_LOG
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
#endif

        BBDSparseMatrix A(M);

#ifdef TIME_LOG
#ifdef _OPENMP
        tper = omp_get_wtime() - t;
        cout << "Время на декомпозицию " << tper << endl;
        t = omp_get_wtime();
#else
        tper = double((clock()-t))/CLOCKS_PER_SEC;
        cout << "Время на декомпозицию " << tper << endl;
        t = clock();
#endif
#endif

        FactorizedBBDSparseMatrix FM;
        MatrixOperations::BlockMatrixFactorization(A, FM);

#ifdef TIME_LOG
#ifdef _OPENMP
        tdec = omp_get_wtime() - t;
        cout << "Время на разложение " << tdec << endl;
        t = omp_get_wtime();
#else
        tdec = double((clock()-t))/CLOCKS_PER_SEC;
        cout << "Время на разложение " << tdec << endl;
        t = clock();
#endif
#endif

        MatrixOperations::Solve(FM, B, X);

#ifdef TIME_LOG
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

#ifdef CSV_LOG
        out << /*tper << "; " << tdec << "; " << tsol << "; " <<*/ tres << "; ";
#endif
#endif

        //X.save2file("solution.txt");

        cout << "Ошибка: " << MatrixOperations::MaxResidual(M, B, X) << endl << endl;

        //cout << X.V[1] << endl;
        //cout << X.V[2] << endl;

    }
    else if (type == 3) //Решение с плотной матрицей
    {

        Matrix M("matrix.txt");

#ifdef TIME_LOG
#ifdef _OPENMP
        omp_set_num_threads(threads);
        double t = omp_get_wtime();
        double t1 = t;
#else
        time_t t = clock();
        time_t t1 = t;
#endif
#endif

        LUPM T;
        MatrixOperations::LUTriang(M, T);

#ifdef TIME_LOG
#ifdef _OPENMP
        cout << "Время на разложение " << omp_get_wtime() - t << endl;
        t = omp_get_wtime();
#else
        cout << "Время на разложение " << double((clock()-t))/CLOCKS_PER_SEC << endl;
        t = clock();
#endif
#endif

        Vector B("vector.txt"), X;
        MatrixOperations::Solve(T, B, X);

#ifdef TIME_LOG
#ifdef _OPENMP
        cout << "Время на решение " << omp_get_wtime() - t << endl;
        cout << "Общее время " << omp_get_wtime() - t1 << endl;
#else
        cout << "Время на решение " << double((clock()-t))/CLOCKS_PER_SEC << endl;
        cout << "Общее время " << double((clock()-t1))/CLOCKS_PER_SEC << endl;
#endif
#endif

        cout << "Ошибка: " << MatrixOperations::MaxResidual(M, B, X) << endl;

        cout << X.V[1] << endl;
        cout << X.V[2] << endl;
    }

    return 0;
}
