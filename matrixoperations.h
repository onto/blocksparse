#ifndef MATRIXOPERATIONS_H
#define MATRIXOPERATIONS_H

#include <cstring>
#include <cmath>
#include <omp.h>

#include "sparsematrix.h"
#include "blockmatrix.h"
#include "matrix.h"
#include "matrixcontainers.h"

#ifdef __NVCC__
#include "cudaoperations.h"
#endif

#define FAST

class MatrixOperations
{
public:
    static bool LUTriang(SparseMatrix &M, LUPS &_M, double pivRel = 0.001);
    static bool LUTriang(Matrix &M, LUPM &_M);

    static bool Solve(LUPS &M, Vector &B, Vector &X);
    static bool Solve(LUPM &M, Vector &B, Vector &X);

    static bool Inverse(SparseMatrix &M, Matrix &_M);
    static bool Inverse(Matrix &M, Matrix &_M);

    static bool BlockMatrixFactorization(BBDSparseMatrix &M,
                                         FactorizedBBDSparseMatrix &_M);
    static bool Solve(FactorizedBBDSparseMatrix &M, Vector &B, Vector &X);

    static double MaxResidual(SparseMatrix &M, Vector &B, Vector &X);
    static double MaxResidual(Matrix &M, Vector &B, Vector &X);
};

bool MatrixOperations::LUTriang(SparseMatrix &M, LUPS &_M, double pivRel)
{
    // Матрицы L и U храним в одной _M.U, но _M.U.F будет означать начало U,
    // а _M.LF будет означать начало L.

    if (M.H != M.W) return false;

    size_t H = M.H;

    size_t j;
    double valabs;

    _M.U = M;
    _M.LF.assign(H+1, SPARSE_END);
    _M.P.resize(H+H+1);
    _M.Pt.resize(H+H+1);
    for (size_t i = 1; i <= H; ++i) _M.P[i] = _M.P[H+i] = i;

    // Массив "указателей" на последний элемент матрицы L, чтобы можно было
    // в неё переносить элементы из U
    std::vector<size_t> LE(H+1, SPARSE_END);

    std::vector<size_t> Rfill(H+1,0), Cfill(H+1,0);

    std::vector<size_t> &F = _M.U.F;
    std::vector<size_t> &LF = _M.LF;
    std::vector<size_t> &C = _M.U.C;
    std::vector<size_t> &N = _M.U.N;
    std::vector<double> &V = _M.U.V;

    //Предрассчитаем количество элементов в строках и стообцах
    for (size_t k = 1; k <= H; ++k)
    {
        for (size_t q = F[k]; q != SPARSE_END; q = N[q])
        {
            ++Rfill[k];
            ++Cfill[C[q]];
        }
    }

    std::vector<double> Norm(H+1);

    for (size_t k = 1; k <= H; ++k)
    {
        //Норма активной подматрицы
        memset(&Norm[0],0,(H+1)*sizeof(double));
        for (size_t i = k; i <= H; ++i)
        {
            for (size_t q = F[i]; q != SPARSE_END; q = N[q])
            {
                if ((valabs = fabs(V[q])) > Norm[(j = C[q])])
                    Norm[j] = valabs;
            }
        }

        size_t optgrowth = H*H, opti = k, optj = k;
        double optv = 0;

        for (size_t i = k; i <= H; ++i)
        {
            for (size_t q = F[i]; q != SPARSE_END; q = N[q])
            {
                j = C[q];
                if ((valabs = fabs(V[q])) >= pivRel * Norm[j])
                {
                    size_t growth = (Rfill[i] - 1) * (Cfill[j] - 1);

                    if ((growth < optgrowth) ||
                            ((growth == optgrowth) && (valabs > optv)))
                    {
                        optgrowth = growth;
                        opti = i;
                        optj = j;
                        optv = valabs;
                        if (optgrowth == 0) break;
                    }
                }                
            }
            if (optgrowth == 0) break;
        }

        // Перестановки строк и столбцов
        // при этом перестановки столбцов U не меняют матрицу L
        size_t t;
        if (opti != k)
        {
            _M.U.swapRow(opti, k);
            t = _M.P[opti]; _M.P[opti] = _M.P[k]; _M.P[k] = t;
            t = LF[opti]; LF[opti] = LF[k]; LF[k] = t;
            t = LE[opti]; LE[opti] = LE[k]; LE[k] = t;
            t = Rfill[opti]; Rfill[opti] = Rfill[k]; Rfill[k] = t;
        }
        if (optj != k)
        {
            _M.U.swapCol(optj, k);
            t = _M.P[H+optj]; _M.P[H+optj] = _M.P[H+k]; _M.P[H+k] = t;
            t = Cfill[optj]; Cfill[optj] = Cfill[k]; Cfill[k] = t;
        }

        std::vector<size_t> kJ;
        std::vector<double> kV;

        size_t q = F[k];
        double diag = V[q]; if (diag == 0.0) V[q] = diag = 1e-20;

        // Преобразуем строку в U
        // и запомним столбцовые индексы и значения в строке
        for (q = N[q]; q != SPARSE_END; q = N[q])
        {
            size_t j = C[q];
            kJ.push_back(j);
            kV.push_back(V[q]/diag);

            --Cfill[j]; //Уменьшаем, т.к. убираем элемент в ведущей строке.
        }
        size_t Nz = kJ.size();

        for (size_t i = k+1; i <= H; ++i)
        {
            j = C[F[i]];
            // Если в строке i нет элемента в веущем столбце, то идем к следующей.
            if (j != k) continue;

            --Rfill[i]; //Уменьшаем, т.к. убираем элемент в ведущем столбце.

            size_t j1, j2;
            size_t p = F[i];
            size_t q = N[p];

            double bdiag = V[p];
            V[p] /= diag;

            for (size_t l = 0; l < Nz; ++l)
            {
                j1 = kJ[l];
                while (q != SPARSE_END)
                {
                    j2 = C[q];
                    if (j2 < j1) { p = q; q = N[q]; continue; }

                    if (j2 == j1)
                    {
                        V[q] -= kV[l]*bdiag;
                    }
                    else
                    {
                        // Тут уже нужно добавить элемент
                        C.push_back(j1);
                        size_t p1 = C.size()-1;
                        N[p] = p1; // предыдущему ссылку на этот
                        N.push_back(q);  // этому ссылку на следующий
                        V.push_back(-kV[l]*bdiag);
                        p = q; q = p1;
                        // больше вроде ничего делать не нужно,
                        // здесь могут попасться только "вторые" элементы в
                        // строке и F менять не нужно.

                        //Увеличим счетчики элементов
                        ++Rfill[i];
                        ++Cfill[j1];
                    }
                    break;
                }

                if (q == SPARSE_END)
                {
                    // Добавим элементы в конец
                    C.push_back(j1);
                    size_t p1 = C.size()-1;
                    N[p] = p1; // предыдущему ссылку на этот
                    N.push_back(q);  // этому ссылку на следующий
                    V.push_back(-kV[l]*bdiag);
                    p = q; q = p1;

                    //Увеличим счетчики элементов
                    ++Rfill[i];
                    ++Cfill[j1];
                }
            }
            // и вот тут элемент в ведущем столбце отходит матрице L
            if (LF[i] == SPARSE_END)
                LF[i] = F[i];
            else
                N[LE[i]] = F[i];
            LE[i] = F[i];
            F[i] = N[F[i]];
            N[LE[i]] = SPARSE_END;
        }
    }

    for (size_t i = 1; i <= H; ++i)
    {
        _M.Pt[_M.P[i]] = i;
        _M.Pt[H+_M.P[H+i]] = i;
    }

    return true;
}

bool MatrixOperations::Solve(LUPS &_M, Vector &B, Vector &X)
{
    size_t H = _M.U.H;

    X.zeros(H);

    // Прямой ход L*y=Pr*B
    for (size_t i = 1; i <= H; ++i) X.V[i] = B.V[_M.P[i]]; //Pr*B;

    for (size_t i = 1; i <= H; ++i)
    {
        for (size_t q = _M.LF[i]; q != SPARSE_END; q = _M.U.N[q])
        {
            X.V[i] -= _M.U.V[q]*X.V[_M.U.C[q]];
        }
    }

    // Обратный ход U*(Pc*x)=y
    for (size_t i = H; i >= 1; --i)
    {
        size_t q = _M.U.F[i];
        double diag = _M.U.V[q];
        for (q = _M.U.N[q]; q != SPARSE_END; q = _M.U.N[q])
        {
            X.V[i] -= _M.U.V[q]*X.V[_M.U.C[q]];
        }
        X.V[i] /= diag;
    }

    std::vector<double> XX = X.V;
    for (size_t i = 1; i <= H; ++i) X.V[i] = XX[_M.Pt[H+i]];

    return true;
}

bool MatrixOperations::Inverse(SparseMatrix &M, Matrix &_M)
{
    LUPS LU;
    if (!LUTriang(M, LU)) return false;

    size_t H = M.H;

    _M.zeros(H, H);

#ifndef __NVCC__

    Vector B(H);
    Vector X;

    // решаем систему LUP * _M = E;
    for (size_t i = 1; i <= H; ++i)
    {
        B.V[i-1] = 0;
        B.V[i] = 1.0;

        if (!Solve(LU, B, X)) return false;

        for (size_t j = 1; j <= H; ++j)
            _M.M[_M.W*(j-1)+i] = X.V[j];
    }

#else

    //Копируем разложенную матрицу на видеокарту
    size_t *cP;     size_t Pbytes = sizeof(size_t)*LU.P.size();
    size_t *cLF;    size_t LFbytes = sizeof(size_t)*LU.LF.size();
    size_t *cC;     size_t Cbytes = sizeof(size_t)*LU.U.C.size();
    size_t *cN;     size_t Nbytes = sizeof(size_t)*LU.U.N.size();
    size_t *cF;     size_t Fbytes = sizeof(size_t)*LU.U.F.size();
    double *cV;     size_t Vbytes = sizeof(double)*LU.U.V.size();

    cudaMalloc((void**)&cP, Pbytes);
    cudaMemcpy(cP, &(LU.P[0]), Pbytes, cudaMemcpyHostToDevice);

    cudaMalloc((void**)&cLF, LFbytes);
    cudaMemcpy(cLF, &(LU.LF[0]), LFbytes, cudaMemcpyHostToDevice);

    cudaMalloc((void**)&cC, Cbytes);
    cudaMemcpy(cC, &(LU.U.C[0]), Cbytes, cudaMemcpyHostToDevice);

    cudaMalloc((void**)&cN, Nbytes);
    cudaMemcpy(cN, &(LU.U.N[0]), Nbytes, cudaMemcpyHostToDevice);

    cudaMalloc((void**)&cF, Fbytes);
    cudaMemcpy(cF, &(LU.U.F[0]), Fbytes, cudaMemcpyHostToDevice);

    cudaMalloc((void**)&cV, Vbytes);
    cudaMemcpy(cV, &(LU.U.V[0]), Vbytes, cudaMemcpyHostToDevice);

    //Выходной вектор с матрицей
    double *cM;     size_t Mbytes = (H*H+1)*sizeof(double);
    cudaMalloc((void**)&cM, Mbytes);
    cudaMemset(cM, 0, Mbytes);

    //Пора запускать ядро
    cudaDeviceProp cProp;
    cudaGetDeviceProperties(&cProp, 0);
    size_t maxThreads = cProp.maxThreadsPerBlock;

    dim3 threads(std::min(H,maxThreads), 1, 1);
    dim3 blocks((H+maxThreads-1)/maxThreads, 1);

    CudaInverse<<<blocks, threads>>>(H, cP, cLF, cC, cN, cF, cV, cM);

    //Эвенты
    cudaEvent_t syncEvent;
    cudaEventCreate(&syncEvent);
    cudaEventRecord(syncEvent, 0);
    cudaEventSynchronize(syncEvent);

    //Заберем матрицу(заодно переставляя строки)
    for (size_t i = 1; i <= H; ++i)
        cudaMemcpy(&(_M.M[H*(i-1)+1]), cM+H*(LU.P[H+i]-1)+1, H*sizeof(double),
                cudaMemcpyDeviceToHost);

    cudaEventDestroy(syncEvent);

    // Уберем за собой
    cudaFree(cP);
    cudaFree(cLF);
    cudaFree(cC);
    cudaFree(cN);
    cudaFree(cF);
    cudaFree(cV);
    cudaFree(cM);

#endif

    return true;
}

bool MatrixOperations::LUTriang(Matrix &M, LUPM &_M)
{
    if (M.H != M.W) return false;

    size_t H = M.H;

    _M.M = M;
    _M.P.resize(H+H+1);
    _M.Pt.resize(H+H+1);
#pragma omp parallel for
    for (size_t i = 1; i <= H; ++i) _M.P[i] = _M.P[H+i] = i;

    for (size_t k = 1; k < H; ++k)
    {
        double norm = 0, value;
        size_t opti = k, optj = k, t;

        for (size_t i = k; i <= H; ++i)
        {
            size_t h = H*(i-1);
            for (size_t j = k; j <= H; ++j)
            {
                value = fabs(_M.M.M[h+j]);
                if (value > norm)
                {
                    norm = value;
                    opti = i;
                    optj = j;
                }
            }
        }

        if (opti != k)
        {
            _M.M.swapRow(opti, k);
            t = _M.P[opti]; _M.P[opti] = _M.P[k]; _M.P[k] = t;
        }
        if (optj != k)
        {
            _M.M.swapCol(optj, k);
            t = _M.P[H+optj]; _M.P[H+optj] = _M.P[H+k]; _M.P[H+k] = t;
        }

        double diag = _M.M.M[_M.M.W*(k-1)+k];

        size_t h1 = H*(k-1);
#pragma omp parallel for
        for (size_t i = k+1; i <= H; ++i)
        {
            size_t h = H*(i-1);
            _M.M.M[h+k] /= diag;
            double nbdiag = _M.M.M[h+k];
            for (size_t j = k+1; j <= H; ++j)
            {
                _M.M.M[h+j] -= _M.M.M[h1+j] * nbdiag;
            }
        }
    }

    for (size_t i = 1; i <= H; ++i)
    {
        _M.Pt[_M.P[i]] = i;
        _M.Pt[H+_M.P[H+i]] = i;
    }

    return true;
}

bool MatrixOperations::Solve(LUPM &_M, Vector &B, Vector &X)
{
    size_t H = _M.M.H;

    X.zeros(H);

    // Прямой ход L*y=Pr*B
    for (size_t i = 1; i <= H; ++i) X.V[i] = B.V[_M.P[i]]; //Pr*B;

    for (size_t i = 1, h = 0; i <= H; ++i, h = H*(i-1))
    {
        for (size_t j = 1; j < i; ++j)
        {
            X.V[i] -= _M.M.M[h+j]*X.V[j];
        }
    }

    // Обратный ход U*(Pc*x)=y
    for (size_t i = H, h = H*(H-1); i >= 1; --i, h = H*(i-1))
    {
        double diag = _M.M.M[h+i];
        for (size_t j = i+1; j <= H; ++j)
        {
            X.V[i] -= _M.M.M[h+j]*X.V[j];
        }
        X.V[i] /= diag;
    }

    std::vector<double> XX = X.V;
    for (size_t i = 1; i <= H; ++i) X.V[i] = XX[_M.Pt[H+i]];

    return true;
}

bool MatrixOperations::Inverse(Matrix &M, Matrix &_M)
{
    LUPM LU;
    if (!LUTriang(M, LU)) return false;

    size_t H = M.H;

    _M.zeros(H, H);

    Vector B(H);
    Vector X;

    for (size_t i = 1; i <= H; ++i)
    {
        B.V[i-1] = 0;
        B.V[i] = 1.0;

        if (!Solve(LU, B, X)) return false;

        for (size_t j = 1; j <= H; ++j)
            _M(j,i) = X.V[j];
    }
    return true;
}

bool MatrixOperations::BlockMatrixFactorization(BBDSparseMatrix &M,
                                                FactorizedBBDSparseMatrix &_M)
{
    size_t Nb = M.Nb; //Количество диагональных блоков

    _M.Nb = M.Nb;
    _M.N = M.N;
    _M.C = M.C;
    _M.BBDS = M.BBDS;

#ifdef TIME_LOG
#ifdef _OPENMP
    double t = omp_get_wtime(), t1;
#endif
#endif

    //Раскладываем A_i
    _M.Alu.resize(Nb);
#pragma omp parallel for
    for (size_t i = 0; i < Nb; ++i)
    {
        LUTriang(M.A[i], _M.Alu[i]);
    }

#ifdef TIME_LOG
#ifdef _OPENMP
    t1 = omp_get_wtime() - t;
    std::cout << "Разложение A_i-x " << t1 << std::endl;
    t = omp_get_wtime();
#endif
#endif

    //Рассчитываем inv(A_i)*B_i = Bh_i
    for (size_t i = 0; i < Nb; ++i)
        _M.Bh.push_back(Matrix(M.B[i]));

#pragma omp parallel for
    for (size_t i = 0; i < Nb; ++i)
    {
        size_t H = _M.Bh[i].H;
        Vector b(H), x(H);
        for (size_t j = 1; j <= M.BBDS.Is[i]; ++j)
        {
            for (size_t k = 1; k <= H; ++k)
                b[k] = _M.Bh[i](k,j);

            Solve(_M.Alu[i], b, x);

            for (size_t k = 1; k <= H; ++k)
                _M.Bh[i](k, j) = x[k];
        }
    }

#ifdef TIME_LOG
#ifdef _OPENMP
    t1 = omp_get_wtime() - t;
    std::cout << "Нахождение Bh_i " << t1 << std::endl;
    t = omp_get_wtime();
#endif
#endif

    // Рассчитаем матрицу H
    Matrix H(M.D);
    for (size_t i = 0; i < Nb; ++i)
    {
        H -= _M.C[i]*_M.Bh[i];
    }

#ifdef TIME_LOG
#ifdef _OPENMP
    t1 = omp_get_wtime() - t;
    std::cout << "Расчет H " << t1 << std::endl;
    t = omp_get_wtime();
#endif
#endif

    // И раскладываем H
    if (!LUTriang(H, _M.H)) return false;

#ifdef TIME_LOG
#ifdef _OPENMP
    t1 = omp_get_wtime() - t;
    std::cout << "Разложение H " << t1 << std::endl;
#endif
#endif

    return true;
}

bool MatrixOperations::Solve(FactorizedBBDSparseMatrix &M, Vector &B,
                             Vector &X)
{
    size_t Nb = M.Nb, N = M.N;

    std::vector< Vector >bh(Nb);
#pragma omp parallel for
    for (size_t ib = 0; ib < Nb; ++ib)
    {
        size_t b0 = M.BBDS.R[ib], bdim = M.BBDS.R[ib+1]-b0;
        Vector b(bdim);
        for (size_t i = 1; i <= bdim; ++i)
            b[i] = B[M.BBDS.Pr[b0+i]];//С учетом перестановки

        Solve(M.Alu[ib], b, bh[ib]);
    }

    Vector g(M.H.M.H); size_t b0 = M.BBDS.R[Nb];
    for (size_t i = 1; i <= M.H.M.H; ++i)
        g[i] = B[M.BBDS.Pr[b0+i]];

    for (size_t i = 0; i < Nb; ++i)
    {
        g -= M.C[i]*bh[i];
    }

    Vector X_q;
    if (!Solve(M.H, g, X_q)) return false;

    X.resize(N);
    Vector Xt(N);

//#pragma omp parallel for
    for (size_t ib = 0; ib < Nb; ++ib)
    {
        size_t bb = M.BBDS.R[ib]+1, bdim = M.BBDS.R[ib+1]-bb+1;

        Vector X_i = bh[ib]-M.Bh[ib]*X_q;

        memcpy(&(Xt.V[bb]),&(X_i.V[1]),bdim*sizeof(double));
    }

    memcpy(&(Xt.V[M.BBDS.R[Nb]+1]),&(X_q.V[1]),M.H.M.H*sizeof(double));

    for (size_t i = 1; i <= N; ++i)
        X[i] = Xt[M.BBDS.Pct[i]];

    return true;
}

double MatrixOperations::MaxResidual(SparseMatrix &M, Vector &B, Vector &X)
{
    Vector R = B-M*X;

    return R.normInf();
}

double MatrixOperations::MaxResidual(Matrix &M, Vector &B, Vector &X)
{
    Vector R = B-M*X;

    return R.normInf();
}

#endif // MATRIXOPERATIONS_H
