#ifndef CUDAKERNELS_CU
#define CUDAKERNELS_CU

#include <cuda.h>
#include <cuda_runtime.h>

__global__ void CudaInverse(size_t H, size_t *P, size_t *LF, size_t *C,
                            size_t *N, size_t *F, double *V, double *M)
{
    int ind = threadIdx.x + blockIdx.x*blockDim.x + 1;

    if ((ind > H) || (ind < 1)) return;

    // Прямой ход L*y=Pr*B
    for (size_t i = 1; i <= H; ++i)
        if (P[i] == ind)
        {
            M[H*(i-1)+ind] = 1;
            break;
        }

    for (size_t i = 1; i <= H; ++i)
    {
        size_t h = H*(i-1)+ind;
        for (size_t q = LF[i]; q != SPARSE_END; q = N[q])
        {
            M[h] -= V[q]*M[H*(C[q]-1)+ind];
        }
    }

    // Обратный ход U*(Pc*x)=y
    for (size_t i = H; i >= 1; --i)
    {
        size_t q = F[i];
        double diag = V[q];
        size_t h = H*(i-1)+ind;
        for (q = N[q]; q != SPARSE_END; q = N[q])
        {
            M[h] -= V[q]*M[H*(C[q]-1)+ind];
        }
        M[h] /= diag;
    }
}

#endif
