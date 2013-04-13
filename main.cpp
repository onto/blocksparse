#include <iostream>

//#include <blockmatrix.h>
//#include <matrix.h>
#include "sparsematrix.h"
#include "matrixoperations.h"
#include <cstdlib>

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

        cout << M.N.size()-1 << endl;

        LUPS T;

        MatrixOperations::LUTriang(M, T);

        Vector X;
        MatrixOperations::Solve(T, B, X);

        //X.print();
    }
    else
    {
        BlockSparseMatrix M1("matrix_b.txt");
        Vector B("vector.txt");

#ifdef _OPENMP
    omp_set_num_threads(threads);
    double t = omp_get_wtime();
#endif
        FactorizedBlockSparseMatrix FM;

        MatrixOperations::BlockMatrixFactorization(M1, FM);
#ifdef _OPENMP
    t = omp_get_wtime() - t;
    cout << "Время на разложение " << t << endl;
    t = omp_get_wtime();
#endif
        Vector X;
        MatrixOperations::Solve(FM, B, X);

#ifdef _OPENMP
    t = omp_get_wtime() - t;
    cout << "Время на решение " << t << endl;
#endif
        //X.print();
    }



//    M.print();
    //cout << '\n';

////    M.swapCol(2,3);

////    M.print();
////    cout << '\n';

//    LUPS<double> T;

//    if (!MatrixOperations<double>::LUTriang(M, T))
//    {
//        cout << "LU :(" << endl;
//    }

//    T.printL();
//        cout << '\n';
//    T.printU();
//        cout << '\n';

//    Matrix<double> Minv;

//    MatrixOperations<double>::Inverse(M,Minv);

//    Minv.print();
//    //cout << '\n';

//    Matrix<double> Minvinv;

//    MatrixOperations<double>::Inverse(Minv, Minvinv);

//    (M-Minvinv).print();
//    cout << '\n';

////    Vector<double> X, B(11);

////    B.set(1,2.9);
////    B.set(2,3.8);

////    MatrixOperations<double>::Solve(T, B, X);

////    X.print();


    //M.A[681].save2fileold("a681old.txt");

    //SparseMatrix<double> M("m6.txt");
//    SparseMatrix<double> M("a681.txt");

//    LUPS<double> T;

//    if (!MatrixOperations<double>::LUTriang(M, T))
//    {
//        cout << "LU :(" << endl;
//    }

//    T.printL();
//        cout << '\n';
//    T.printU();
//        cout << '\n';


    return 0;
}
