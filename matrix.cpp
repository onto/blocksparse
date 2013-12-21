#include "matrix.h"

Matrix::Matrix()
{
    H = 0;
    W = 0;
}
Matrix::Matrix(const Matrix &_M)
{
    H = _M.H;
    W = _M.W;
    M = _M.M;
}
Matrix::Matrix(size_t h, size_t w)
{
    zeros(h,w);
}
Matrix::Matrix(size_t h)
{
    zeros(h,h);
}
Matrix::Matrix(const char *file)
{
    std::ifstream in(file);

    in >> H >> W;

    zeros(H, W);

    double x;
    size_t c;

    for (size_t r = 1; r <= H; ++r)
    {
        while (true)
        {
            in >> c;

            if (c != SPARSE_END) {
                in >> x;
                set(r, c, x);
            } else {
                break;
            }
        }
    }

    in.close();
}

Matrix::Matrix(const SparseMatrix &S)
{
    zeros(S.H, S.W);

    for (size_t i = 1; i <= H; ++i)
    {
        size_t h = W*(i-1);
        for (size_t q = S.F[i]; q != SPARSE_END; q = S.N[q])
        {
            M[h+S.C[q]] = S.V[q];
        }
    }
}

double Matrix::get(size_t row, size_t col) const
{
    return M[W*(row-1)+col];
}

void Matrix::set(size_t row, size_t col, double value)
{
    M[W*(row-1)+col] = value;
}

void Matrix::zeros(size_t h, size_t w)
{
    H = h;
    W = w;
    M.resize(H*W+1);
    memset(&M[0], 0, (H*W+1)*sizeof(double));
}

void Matrix::print() const
{
    for (size_t i = 1; i <= H; ++i)
    {
        for (size_t j = 1; j <= W; ++j)
        {
            if (fabs(get(i,j)) < 1e-10)
                std::cout << "0\t";
            else
                std::cout << get(i,j) << '\t';
        }
        std::cout << ";\n";
    }
}

void Matrix::swapRow(size_t r1, size_t r2)
{
    double *t = new double[W];

    memcpy(t, &M[W*(r1-1)+1], W*sizeof(double));
    memcpy(&M[W*(r1-1)+1], &M[W*(r2-1)+1], W*sizeof(double));
    memcpy(&M[W*(r2-1)+1], t, W*sizeof(double));

    delete[] t;
}

void Matrix::swapCol(size_t c1, size_t c2)
{
#pragma omp parallel for
    for (size_t i = 1; i <= H; ++i)
    {
        size_t h1 = W*(i-1)+c1, h2 = W*(i-1)+c2;
        double t = M[h1];
        M[h1] = M[h2];
        M[h2] = t;
    }
}

Matrix operator +(const Matrix &M1, const Matrix &M2)
{
    if ((M1.H != M2.H) || (M1.W != M2.W)) {} // надо бы бросить исключение

    Matrix OM(M1);

#pragma omp parallel for
    for (size_t i = 1; i <= OM.M.size(); ++i)
        OM.M[i] += M2.M[i];

    return OM;
}

Matrix operator -(const Matrix &M1, const Matrix &M2)
{
    if ((M1.H != M2.H) || (M1.W != M2.W)) {} // надо бы бросить исключение

    Matrix OM(M1);

#pragma omp parallel for
    for (size_t i = 1; i <= OM.M.size(); ++i)
        OM.M[i] -= M2.M[i];

    return OM;
}

Matrix& operator -=(Matrix &M1, const Matrix &M2)
{
    if ((M1.H != M2.H) || (M1.W != M2.W)) {} // надо бы бросить исключение

#pragma omp parallel for
    for (size_t i = 1; i < M1.M.size(); ++i)
        M1.M[i] -= M2.M[i];

    return M1;
}

Matrix& operator +=(Matrix &M1, const Matrix &M2)
{
    if ((M1.H != M2.H) || (M1.W != M2.W)) {} // надо бы бросить исключение

#pragma omp parallel for
    for (size_t i = 1; i < M1.M.size(); ++i)
        M1.M[i] += M2.M[i];

    return M1;
}

Matrix operator *(const Matrix &M1, const Matrix &M2)
{
    if (M1.W != M2.H) {} // надо бы бросить исключение

    Matrix OM(M1.H, M2.W);

#pragma omp parallel for
    for (size_t i = 1; i <= M1.H; ++i)
    {
        size_t h1 = M1.W*(i-1);
        for (size_t j = 1; j <= M2.W; ++j)
        {
            size_t h2 = OM.W*(i-1)+j;
            for (size_t k = 1; k <= M1.W; ++k)
                OM.M[h2] += M1.M[h1+k] * M2.M[M2.W*(k-1)+j];
        }
    }

    return OM;
}

Vector operator *(const Matrix &M, const Vector &V)
{
    if (M.W != V.H) {} // надо бы бросить исключение

    Vector OV(M.H);

#pragma omp parallel for
    for (size_t i = 1; i <= M.H; ++i)
    {
        size_t h = M.W*(i-1);
        for (size_t j = 1; j <= M.W; ++j)
            OV.V[i] += M.M[h+j] * V.V[j];
    }

    return OV;
}

Matrix operator *(const SparseMatrix &S, const Matrix &M)
{
    if (S.W != M.H) {} // надо бы бросить исключение

    Matrix OM(S.H,M.W);

#pragma omp parallel for
    for (size_t i = 1; i <= OM.H; ++i)
    {
        size_t h = OM.W*(i-1);
        for (size_t q = S.F[i]; q != SPARSE_END; q = S.N[q])
        {
            size_t h1 = M.W*(S.C[q]-1);
            double v = S.V[q];
            for (size_t k = 1; k <= OM.W; ++k)
            {
                //OM[i,k] += S[i,j]*M[j,k]
                OM.M[h+k] += v * M.M[h1+k];
            }
        }
    }

    return OM;
}

Matrix operator *(const Matrix &M, const SparseMatrix &S)
{
    if (M.W != S.H) {} // надо бы бросить исключение

    Matrix OM(M.H,S.W);

    for (size_t i = 1; i <= M.W; ++i)
    {
        for (size_t q = S.F[i]; q != SPARSE_END; q = S.N[q])
        {
            size_t j = S.C[q];
            double v = S.V[q];
#pragma omp parallel for
            for (size_t k = 1; k <= OM.H; ++k)
            {
                //OM[k,j] += M[k,i]*S[i,j]
                OM.M[OM.W*(k-1)+j] += M.M[M.W*(k-1)+i] * v;
            }
        }
    }

    return OM;
}

Vector operator *(const SparseMatrix &S, const Vector &V)
{
    if (S.W != V.H) {} // надо бы бросить исключение

    Vector OV(S.H);

#pragma omp parallel for
    for (size_t i = 1; i <= S.H; ++i)
        for (size_t q = S.F[i]; q != SPARSE_END; q = S.N[q])
            OV.V[i] += S.V[q] * V.V[S.C[q]];

    return OV;
}

Matrix operator +(const SparseMatrix &S, const Matrix &M)
{
    if ((S.H != M.H) || (S.W != M.W)) {} // надо бы бросить исключение

    Matrix OM(M);

#pragma omp parallel for
    for (size_t i = 1; i <= S.H; ++i)
    {
        size_t h = OM.W*(i-1);
        for (size_t q = S.F[i]; q != SPARSE_END; q = S.N[q])
            OM.M[h+S.C[q]] += S.V[q];
    }

    return OM;
}

Matrix operator +(const Matrix &M, const SparseMatrix &S)
{
    return S + M;
}

Matrix& operator +=(Matrix &M, const SparseMatrix &S)
{
    if ((S.H != M.H) || (S.W != M.W)) {} // надо бы бросить исключение

#pragma omp parallel for
    for (size_t i = 1; i <= S.H; ++i)
    {
        size_t h = M.W*(i-1);
        for (size_t q = S.F[i]; q != SPARSE_END; q = S.N[q])
            M.M[h+S.C[q]] += S.V[q];
    }

    return M;
}

Matrix operator -(const SparseMatrix &S, const Matrix &M)
{
    if ((S.H != M.H) || (S.W != M.W)) {} // надо бы бросить исключение

    Matrix OM(M);

#pragma omp parallel for
    for (size_t i = 1; i < OM.M.size(); ++i) OM.M[i] = -OM.M[i];

#pragma omp parallel for
    for (size_t i = 1; i <= S.H; ++i)
    {
        size_t h = OM.W*(i-1);
        for (size_t q = S.F[i]; q != SPARSE_END; q = S.N[q])
            OM.M[h+S.C[q]] += S.V[q];
    }

    return OM;
}

Matrix operator -(const Matrix &M, const SparseMatrix &S)
{
    if ((S.H != M.H) || (S.W != M.W)) {} // надо бы бросить исключение

    Matrix OM(M);

#pragma omp parallel for
    for (size_t i = 1; i <= S.H; ++i)
    {
        size_t h = OM.W*(i-1);
        for (size_t q = S.F[i]; q != SPARSE_END; q = S.N[q])
            OM.M[h+S.C[q]] -= S.V[q];
    }

    return OM;
}

Vector::Vector()
{
    H = 0;
}
Vector::Vector(const Vector &V1)
{
    V = V1.V;
    H = V1.H;
}
Vector::Vector(size_t h)
{
    H = h;
    V.resize(H+1);
}

Vector::Vector(std::vector<double> &V1)
{
    H = V1.size()-1;
    V = V1;
}

Vector::Vector(const char *file)
{
    std::ifstream in(file);

    V.push_back(0);

    double x;
    while (in >> x)
        V.push_back(x);

    H = V.size()-1;
}

inline double Vector::get(size_t i) const
{
    return V[i];
}

inline void Vector::set(size_t i, double value)
{
    V[i] = value;
}

void Vector::zeros(size_t h)
{
    H = h;
    V.resize(H+1);
    memset(&V[0], 0, (H+1)*sizeof(double));
}

void Vector::resize(size_t h)
{
    zeros(h);
}

void Vector::print() const
{
    for (size_t i = 1; i <= H; ++i)
        std::cout << V[i] << '\n';
}

void Vector::permute(std::vector<size_t>& P)
{
    std::vector<double> Vtemp(H+1);
    for (size_t i = 1; i <= H; ++i)
        Vtemp[P[i]] = V[i];

    V = Vtemp;
}

Vector operator +(const Vector &V1, const Vector &V2)
{
    if (V1.H != V2.H) {}

    Vector OV(V1.H);

#pragma omp parallel for
    for (size_t i = 0; i <= V1.H; ++i)
        OV.V[i] = V1.V[i] + V2.V[i];

    return OV;
}

Vector operator -(const Vector &V1, const Vector &V2)
{
    if (V1.H != V2.H) {}

    Vector OV(V1.H);

#pragma omp parallel for
    for (size_t i = 0; i <= V1.H; ++i)
        OV.V[i] = V1.V[i] - V2.V[i];

    return OV;
}

Vector& operator -=(Vector &V1, const Vector &V2)
{
    if (V1.H != V2.H) {}

#pragma omp parallel for
    for (size_t i = 0; i <= V1.H; ++i)
        V1.V[i] -= V2.V[i];

    return V1;
}

Vector& operator +=(Vector &V1, const Vector &V2)
{
    if (V1.H != V2.H) {}

#pragma omp parallel for
    for (size_t i = 0; i <= V1.H; ++i)
        V1.V[i] += V2.V[i];

    return V1;
}


