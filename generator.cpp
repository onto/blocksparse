#include <iostream>
#include <fstream>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <vector>

#include "sparsematrix.h"

#define RANDN ((rand()%2 == 1)?1.0:-1.0)*(1.0+rand())/(1.0+rand())

using namespace std;

int main(int argc, char *argv[4])
{
    if (argc == 2) //Генерация разреженной симметричной матрицы
    {
        size_t N = (size_t)atoi(argv[1]);

        SparseMatrix S(N, N);

        size_t p = N / 4 + 1;

        srand(time(NULL));

        size_t q;

        for (size_t i = 1; i <= N; ++i)
        {
            q = (i==N)?N:rand()%(N-i) + i+1;
            for (size_t j = i+1; j <= N; ++j)
            {
                if ((rand()%p == 0) || (j == q))
                {
                    S.set(i, j, RANDN);
                    S.set(j, i, RANDN);
                }
            }
            //if (rand()%p == 0)
            //{
                S.set(i, i, RANDN);
            //}
        }

        //S.print();

        S.save2file("matrix.txt");

        ofstream V("vector.txt");
        for (size_t i = 0; i < N; ++i)
            V << RANDN << '\t';
        V << endl;

        //cout << "gen OK" << endl;

        return 0;
    }

    if (argc < 4)
    {
        cout << "Need dimension, block maxsize, minsize" << endl;
        return 0;
    }

    size_t N = (size_t)atoi(argv[1]);

    srand(time(NULL));

    ofstream M1("matrix_b.txt");
    ofstream M2("matrix_s.txt");
    ofstream V("vector.txt");

    int Max = atoi(argv[2]);
    int Min = atoi(argv[3]);

    vector<size_t> R;

    size_t q = 0, t;
    int Q;
    while(true)
    {
        t = rand()%(Max-Min+1)+Min;
        q += t;
        if (q > N) {t -= (q-N); q = N;}
        if (N-q < 2) {t += (N-q); q = N;}
        R.push_back(t);
        M1 << t << '\t';
        if (q == N) break;
    }
    M1 << 0 << endl;
    M2 << N << " " << N << endl;

    q = 1;
    for (size_t k = 0; k < R.size()-1; ++k)
    {
        for (size_t i = q; i < q+R[k]; ++i)
        {
            Q = int((R[k])/4);
            if (Q == 0) Q = 1;

            //Строка A
            for (size_t j = q; j < q+R[k]; ++j)
            {
                double num = RANDN;
                if (i == j)
                {
                    M1 << j << " " << num << '\t';
                    M2 << j << " " << num << '\t';
                }
                else
                {
                    int c = rand()%Q;
                    if (c == 0)
                    {
                        M1 << j << " " << num << '\t';
                        M2 << j << " " << num << '\t';
                    }
                }
            }

            Q = int((R[k]+R.back())/4);
            if (Q == 0) Q = 1;
            //Строка B
            for (size_t j = N-R.back()+1; j <= N; ++j)
            {
                double num = RANDN;
                int c = rand()%Q;
                if (c == 0)
                {
                    M1 << j << " " << num << '\t';
                    M2 << j << " " << num << '\t';
                }
            }
            M1 << 0 << endl;
            M2 << 0 << endl;
        }
        q += R[k];
    }

    //Строки С и Q
    for (size_t i = 0; i < R.back(); ++i)
    {
        q = 1;
        for (size_t k = 0; k < R.size()-1; ++k)
        {
            Q = int((R[k]+R.back())/4);
            if (Q == 0) Q = 1;
            for (size_t j = q; j < q+R[k]; ++j)
            {
                double num = RANDN;
                int c = rand()%Q;
                if (c == 0)
                {
                    M1 << j << " " << num << '\t';
                    M2 << j << " " << num << '\t';
                }
            }

            q += R[k];
        }

        Q = int((R.back())/4);
        if (Q == 0) Q = 1;
        for (size_t j = N-R.back()+1, k = 0; j <= N; ++j, ++k)
        {
            double num = RANDN;
            if (i == k)
            {
                M1 << j << " " << num << '\t';
                M2 << j << " " << num << '\t';
            }
            else
            {
                int c = rand()%Q;
                if (c == 0)
                {
                    M1 << j << " " << num << '\t';
                    M2 << j << " " << num << '\t';
                }
            }
        }
        M1 << 0 << endl;
        M2 << 0 << endl;
    }

    for (size_t i = 0; i < N; ++i)
        V << RANDN << '\t';
    V << endl;

    return 0;
}
