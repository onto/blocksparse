#ifndef MATRIXOPERATIONS_H
#define MATRIXOPERATIONS_H

#include <cstring>
#include <cmath>
#include <omp.h>

#include <sparsematrix.h>
#include <blockmatrix.h>
#include <matrix.h>


template <typename T> struct LUPS;
template <typename T> struct LUPM;

template <typename T>
struct LUPS
{
    SparseMatrix<T> U;
    std::vector<size_t> LF;
    std::vector<size_t> P;

    void printL() const
    {
        for (size_t i = 1; i <= U.H; ++i)
        {
            size_t j = 1, jj;
            for (size_t q = LF[i]; q != SPARSE_END; q = U.N[q])
            {
                jj = U.C[q];
                while (j < jj) { std::cout << "0" << '\t'; ++j; }
                std::cout << U.V[q] << '\t'; ++j;
            }
            while (j <= U.W) { std::cout << "0" << '\t'; ++j; }
            std::cout << ";\n";
        }
    }

    void printU() const
    {
        U.print();
    }
};

template <typename T>
struct LUPM
{
    Matrix<T> M;
    std::vector<size_t> P;

    void printL() const
    {
        for (size_t i = 1; i <= M.H; ++i)
        {
            for (size_t j = 1; j < i; ++j)
                std::cout << M.get(i,j) << '\t';
            for (size_t j = i; j <= M.H; ++j)
                std::cout << "*\t";
            std::cout << '\n';
        }
    }

    void printU() const
    {
        for (size_t i = 1; i <= M.H; ++i)
        {
            for (size_t j = 1; j < i; ++j)
                std::cout << "*\t";
            for (size_t j = i; j <= M.H; ++j)
                std::cout << M.get(i,j) << '\t';
            std::cout << '\n';
        }
    }
};

template <typename T>
class MatrixOperations
{
public:
    static bool LUTriang(SparseMatrix<T> &M, LUPS<T> &_M, double pivRel = 0.001);
    static bool LUTriang(Matrix<T> &M, LUPM<T> &_M);

    static bool Solve(LUPS<T> &M, Vector<T> &B, Vector<T> &X);
    static bool Solve(LUPM<T> &M, Vector<T> &B, Vector<T> &X);

    static bool Inverse(SparseMatrix<T> &M, Matrix<T> &_M);
    static bool Inverse(Matrix<T> &M, Matrix<T> &_M);

    static bool BlockMatrixFactorization(BlockSparseMatrix<T> &M, FactorizedBlockSparseMatrix<T> &_M);
    static bool Solve(FactorizedBlockSparseMatrix<T> &M, Vector<T> &B, Vector<T> &X);
};

template <typename T>
bool MatrixOperations<T>::LUTriang(SparseMatrix<T> &M, LUPS<T> &_M, double pivRel)
{
    // Матрицы L и U храним в одной _M.U, но _M.U.F будет означать начало U,
    // а _M.LF будет означать начало L.

    if (M.H != M.W) return false;

    size_t H = M.H;

    _M.U = M;
    _M.LF.assign(H+1, SPARSE_END);
    _M.P.resize(H+H+1);
    for (size_t i = 1; i <= H; ++i) _M.P[i] = i;
    for (size_t i = 1; i <= H; ++i) _M.P[H+i] = i;

    // Массив "указателей" на последний элемент матрицы L, чтобы можно было
    // в неё переносить элементы из U
    std::vector<size_t> LE(H+1, SPARSE_END);

//    size_t *N = &(_M.U.N[0]), *C = &(_M.U.C[0]), *F = &(_M.U.F[0]),
//            *P = &(_M.P[0]), *LF = &(_M.LF[0]);
//    T *V = &(_M.U.V[0]);

    std::vector<size_t> Rfill, Cfill;

    for (size_t k = 1; k <= H; ++k)
    {
        Rfill.assign(H+1, 0);
        Cfill.assign(H+1, 0);

        double norm = 0;
        for (size_t i = k; i <= H; ++i)
        {
            for (size_t q = _M.U.F[i]; q != SPARSE_END; q = _M.U.N[q])
            {
                double t = fabs(_M.U.V[q]);
                if (t > norm)
                    norm = t;

                ++Rfill[i];
                ++Cfill[_M.U.C[q]];
            }
        }

        size_t optgrowth = H*H, opti = k, optj = k;
        double optv = 0;

        for (size_t i = k; i <= H; ++i)
        {
            //if ((Rfill[i] == 0) || (Cfill[i])) return false;

            for (size_t q = _M.U.F[i]; q != SPARSE_END; q = _M.U.N[q])
            {
                double val = fabs(_M.U.V[q]);
                if (val >= pivRel * norm)
                {

                    size_t j = _M.U.C[q];
                    size_t growth = (Rfill[i] - 1) * (Cfill[j] - 1);

                    if ((growth < optgrowth) ||
                            ((growth == optgrowth) && (val > optv)))
                    {
                        optgrowth = growth;
                        opti = i;
                        optj = j;
                        optv = val;
                    }
                }
            }
        }

//        std::cout << "BEFORE SWAP\n";
//        _M.printL();
//        std::cout << '\n';
//        _M.printU();
//        std::cout << '\n';

        // Перестановки строк и столбцов
        // при этом перестановки столбцов U не меняют матрицу L
        size_t t;
        if (opti != k)
        {
            _M.U.swapRow(opti, k);
            t = _M.P[opti]; _M.P[opti] = _M.P[k]; _M.P[k] = t;
            t = _M.LF[opti]; _M.LF[opti] = _M.LF[k]; _M.LF[k] = t;
            t = LE[opti]; LE[opti] = LE[k]; LE[k] = t;
        }
        if (optj != k)
        {
            _M.U.swapCol(optj, k);
            t = _M.P[H+optj]; _M.P[H+optj] = _M.P[H+k]; _M.P[H+k] = t;
        }

//        std::cout << "SWAP\n";
//        _M.printL();
//        std::cout << '\n';
//        _M.printU();
//        std::cout << '\n';

        std::vector<size_t> kJ;
        std::vector<T> kV;

        T diag = _M.U.V[_M.U.F[k]];

        // Преобразуем строку в U
        // и запомним столбцовые индексы и значения в строке
        for (size_t q = _M.U.N[_M.U.F[k]]; q != SPARSE_END; q = _M.U.N[q])
        {
            kJ.push_back(_M.U.C[q]);
            kV.push_back(_M.U.V[q]/diag);
        }

        for (size_t i = k+1; i <= H; ++i)
        {
            size_t j = _M.U.C[_M.U.F[i]];
            // Если в строке i нет диагонального, то идем к следующей.
            if (j != k) continue;

            size_t j1, j2;
            size_t p = _M.U.F[i];
            size_t q = _M.U.N[p];

            T bdiag = _M.U.V[p];

            _M.U.V[p] /= diag;

            for (size_t l = 0; l < kJ.size(); ++l)
            {
                j1 = kJ[l];
                while (q != SPARSE_END)
                {

                    j2 = _M.U.C[q];
                    if (j2 < j1) { p = q; q = _M.U.N[q]; continue; }

                    if (j2 == j1)
                    {
                        _M.U.V[q] -= kV[l]*bdiag;
                    }
                    else
                    {
                        // Тут уже нужно добавить элемент
                        _M.U.C.push_back(j1);
                        _M.U.N[p] = _M.U.C.size()-1; // предыдущему ссылку на этот
                        _M.U.N.push_back(q);  // этому ссылку на следующий
                        _M.U.V.push_back(-kV[l]*bdiag);
                        p = q; q = _M.U.C.size()-1;
                        // больше вроде ничего делать не нужно,
                        // здесь могут попасться только "вторые" элементы в
                        // строке и F менять не нужно.
                    }
                    break;
                }

                if (q == SPARSE_END)
                {
                    // Добавим элементы в конец
                    _M.U.C.push_back(j1);
                    _M.U.N[p] = _M.U.C.size()-1; // предыдущему ссылку на этот
                    _M.U.N.push_back(q);  // этому ссылку на следующий
                    _M.U.V.push_back(-kV[l]*bdiag);
                    p = q; q = _M.U.C.size()-1;
                }
            }
            // и вот тут элемент в ведущем столбце отходит матрице L
            if (_M.LF[i] == SPARSE_END)
                _M.LF[i] = _M.U.F[i];
            else
                _M.U.N[LE[i]] = _M.U.F[i];
            LE[i] = _M.U.F[i];
            _M.U.F[i] = _M.U.N[_M.U.F[i]];
            _M.U.N[LE[i]] = SPARSE_END;
        }
    }

    return true;
}

template <typename T>
bool MatrixOperations<T>::Solve(LUPS<T> &_M, Vector<T> &B, Vector<T> &X)
{
    size_t H = _M.U.H;

    // Транспонируем P для столбцов
    std::vector<size_t> P(H+1);
    for (size_t i = 1; i <= H; ++i) P[_M.P[H+i]] = i;

    X.H = H;
    X.V.assign(H+1, 0);

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
        T diag = _M.U.V[_M.U.F[i]];
        for (size_t q = _M.U.N[_M.U.F[i]]; q != SPARSE_END; q = _M.U.N[q])
        {
            X.V[i] -= _M.U.V[q]*X.V[_M.U.C[q]];
        }
        X.V[i] /= diag;
    }

    std::vector<double> XX = X.V;
    for (size_t i = 1; i <= H; ++i) X.V[i] = XX[P[i]];

    return true;
}

template <typename T>
bool MatrixOperations<T>::Inverse(SparseMatrix<T> &M, Matrix<T> &_M)
{
    LUPS<T> LU;
    if (!LUTriang(M, LU)) return false;

    size_t H = M.H;

    _M = Matrix<T>(H);

    Vector<T> B(H);
    Vector<T> X;

    // решаем систему LUP * _M = E;
    for (size_t i = 1; i <= H; ++i)
    {

        B.V[i-1] = 0;
        B.V[i] = 1.0;

        if (!Solve(LU, B, X)) return false;

        for (size_t j = 1; j <= H; ++j)
            _M.set(j, i, X.V[j]);
        //memcpy(&(_M.M[1+H*(i-1)]),&(X[1]),sizeof(T)*H);
    }
    return true;
}

template <typename T>
bool MatrixOperations<T>::LUTriang(Matrix<T> &M, LUPM<T> &_M)
{
    if (M.H != M.W) return false;

    size_t H = M.H;

    _M.M = M;
    _M.P.resize(H+H+1);
#pragma omp parallel for
    for (size_t i = 1; i <= H; ++i) _M.P[i] = _M.P[H+i] = i;

    for (size_t k = 1; k < H; ++k)
    {

        double norm = 0, value;
        size_t opti = k, optj = k, t;

        for (size_t i = k; i <= H; ++i)
        {
            for (size_t j = k; j <= H; ++j)
            {
                value = fabs(_M.M.M[H*(i-1)+j]);
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

        T diag = _M.M.get(k,k);

#pragma omp parallel for
        for (size_t i = k+1; i <= H; ++i)
        {
            _M.M.M[H*(i-1)+k] /= diag;
            T nbdiag = _M.M.M[H*(i-1)+k];
            for (size_t j = k+1; j <= H; ++j)
            {
                _M.M.M[H*(i-1)+j] -= _M.M.M[H*(k-1)+j] * nbdiag;
            }
        }

    }

    return true;
}

template <typename T>
bool MatrixOperations<T>::Solve(LUPM<T> &_M, Vector<T> &B, Vector<T> &X)
{
    size_t H = _M.M.H;

    // Транспонируем P для столбцов
    std::vector<size_t> P(H+1);
    for (size_t i = 1; i <= H; ++i) P[_M.P[H+i]] = i;

    X.H = H;
    X.V.assign(H+1, 0);

    // Прямой ход L*y=Pr*B
    for (size_t i = 1; i <= H; ++i) X.V[i] = B.V[_M.P[i]]; //Pr*B;

    for (size_t i = 1; i <= H; ++i)
    {
        for (size_t j = 1; j < i; ++j)
        {
            X.V[i] -= _M.M.M[H*(i-1)+j]*X.V[j];
        }
    }

    // Обратный ход U*(Pc*x)=y
    for (size_t i = H; i >= 1; --i)
    {
        T diag = _M.M.M[H*(i-1)+i];
        for (size_t j = i+1; j <= H; ++j)
        {
            X.V[i] -= _M.M.M[H*(i-1)+j]*X.V[j];
        }
        X.V[i] /= diag;
    }

    std::vector<T> XX = X.V;
    for (size_t i = 1; i <= H; ++i) X.V[i] = XX[P[i]];

    return true;
}

template <typename T>
bool MatrixOperations<T>::Inverse(Matrix<T> &M, Matrix<T> &_M)
{
    LUPM<T> LU;
    if (!LUTriang(M, LU)) return false;

    size_t H = M.H;

    _M = Matrix<T>(H);

    Vector<T> B(H);
    Vector<T> X;

    // решаем систему LUP * _M = E;
    for (size_t i = 1; i <= H; ++i)
    {

        B.V[i-1] = 0;
        B.V[i] = 1.0;

        if (!Solve(LU, B, X)) return false;

        for (size_t j = 1; j <= H; ++j)
            _M.set(j, i, X.V[j]);
        //memcpy(&(_M.M[1+H*(i-1)]),&(X[1]),sizeof(T)*H);
    }
    return true;
}

template <typename T>
bool MatrixOperations<T>::BlockMatrixFactorization(BlockSparseMatrix<T> &M, FactorizedBlockSparseMatrix<T> &_M)
{

    _M.B = M.B;
    _M.C = M.C;
    _M.Q = M.Q;
    _M.R = M.R;

    size_t N = _M.R.size()-1;

    // Инвертируем A_i
    _M.Ainv.resize(N);
    #pragma omp parallel for
    for (size_t i = 0; i < N; ++i)
    {
        Inverse(M.A[i], _M.Ainv[i]);
    }

    // Считаем матрицу H
    Matrix<T> H(_M.R.back());
    for (size_t i = 0; i < N; ++i)
    {
        H -= _M.C[i]*_M.Ainv[i]*_M.B[i];
    }
    H += _M.Q;

    // И раскладываем H
    if (!LUTriang(H, _M.H)) return false;

    return true;
}

template <typename T>
bool MatrixOperations<T>::Solve(FactorizedBlockSparseMatrix<T> &M, Vector<T> &B, Vector<T> &X)
{

    std::vector< Vector<T> >b;

    size_t N = M.R.size()-1;

    size_t q = 1;
    std::vector<size_t> qI(N+2);
    for (size_t i = 0; i <= N; ++i)
    {
        qI[i] = q;
        q += M.R[i];
    }
    qI[N+1] = q;

    b.resize(N+1);
    #pragma omp parallel for
    for (size_t i = 0; i <= N; ++i)
    {
        std::vector<T> t;
        t.assign(B.V.begin()+qI[i], B.V.begin()+qI[i+1]);
        t.insert(t.begin(), 0);
        b[i].H = M.R[i];
        b[i].V = t;
    }

    X = Vector<T>(B.H);

    Vector<T> v = b.back();
    for (size_t i = 0; i < N; ++i)
    {
        v -= M.C[i]*(M.Ainv[i]*b[i]);
    }

    Vector<T> X_q;
    if (!Solve(M.H, v, X_q)) return false;

    for (size_t i = 0; i < N; ++i)
    {
        Vector<T> X_i = M.Ainv[i]*(b[i]-M.B[i]*X_q);
        memcpy(&(X.V[qI[i]]),&(X_i.V[1]),M.R[i]*sizeof(T));
    }

    q = qI[N];
    memcpy(&(X.V[q]),&(X_q.V[1]),M.R.back()*sizeof(T));

    return true;
}




#endif // MATRIXOPERATIONS_H
