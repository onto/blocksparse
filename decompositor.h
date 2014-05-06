#ifndef DECOMPOSITOR_H
#define DECOMPOSITOR_H

#include <cstdlib>
#include <ctime>
#include <vector>
#include <algorithm>
#include <queue>
#include <cstring>
#include <climits>

#include "sparsematrix.h"
#include "blockmatrix.h"

enum DecompositionMethod
{
    SV = 0,
    SVI
};

struct BBDStruct
{
    VectSizet Pr, Pc, Prt, Pct; //Перестановки
    VectSizet R;  //Размеры блоков
    VectSizet Is; //Массив количества элементов в каждом разделителе

    void swap(VectSizet _Pr, VectSizet _Pc, VectSizet _Prt, VectSizet _Pct,
              VectSizet _R, VectSizet _Is)
    {
        Pr.swap(_Pr);
        Pc.swap(_Pc);
        Prt.swap(_Prt);
        Pct.swap(_Pct);
        R.swap(_R);
        Is.swap(_Is);
    }
};

class MatrixBBDPreordering
{
public:
    static void BBDDecompose(SparseMatrix &A, BBDStruct &BBDS, DecompositionMethod method, size_t Nmin, size_t Nmax, size_t Dmax);

private:
    static void MaximumTransversal(SparseMatrix &A, VectSizet &Pr);

    static void SangiovannyVincentelli(VectVectSizet &H, VectVectSizet &D, VectSizet &Is, size_t Nmin, size_t Nmax);
    static void SangiovannyVincentelliImproved(VectVectSizet &H, VectVectSizet &D, VectSizet &Is, size_t Nmin, size_t Nmax, size_t Dmax);
};

void MatrixBBDPreordering::BBDDecompose(SparseMatrix &A, BBDStruct &BBDS, DecompositionMethod method, size_t Nmin, size_t Nmax, size_t Dmax)
{
    size_t N = A.H, N1 = N+1;

#ifdef TIME_LOG
#ifdef _OPENMP
    double t = omp_get_wtime(), t1;
#endif
#endif

    //Для начала избавимся от нулей на диагонали,
    //чтобы диагональные блоки не получались вырожденными
    VectSizet PrMT, PrtMT(N1);
    MaximumTransversal(A, PrMT);

#ifdef TIME_LOG
#ifdef _OPENMP
    t1 = omp_get_wtime() - t;
    std::cout << "MaximumTransversal " << t1 << std::endl;
    t = omp_get_wtime();
#endif
#endif


    for (size_t i = 1; i <= N; ++i) PrtMT[PrMT[i]] = i;

    VectVectSizet H(N1); //Массив хранения структуры матрицы A+A^T

    //Заполняем структуру матрицы A+A^T, применяя к ней перестановку полученную
    //в методе MaximumTransversal
    size_t j;
    for (size_t i = 1; i <= N; ++i)
    {
        size_t ii = PrtMT[i];
        for (size_t q = A.F[i]; q != SPARSE_END; q = A.N[q])
        {
            j = A.C[q];
            if (ii != j)
            {
                H[ii].push_back(j);
                H[j].push_back(ii);
            }
        }
    }

    //Убираем повторяющиеся
    for (size_t i = 1; i <= N; ++i)
    {
        std::sort(H[i].begin(), H[i].end());
        H[i].resize(std::distance(H[i].begin(), std::unique(H[i].begin(), H[i].end())));
    }

    VectVectSizet D;
    VectSizet Is;

    //Рассчитываем перестановку
    switch (method) {
    default:
    case SVI:
        SangiovannyVincentelliImproved(H, D, Is, Nmin, Nmax, Dmax);
        break;
    case SV:
        SangiovannyVincentelli(H, D, Is, Nmin, Nmax);
        break;
    }

    VectSizet Pr, Pc, Prt, Pct;
    VectSizet R;

    //Заполняем массив перестановки столбцов и массив размеров блоков
    Pc.push_back(0);
    R.push_back(0);
    for (size_t i = 1, i_end = D.size(); i < i_end; ++i)
    {
        Pc.insert(Pc.end(), D[i].begin(), D[i].end());
        R.push_back(R.back()+D[i].size());
    }

    if (D[0].size() != 0)
    {
        Pc.insert(Pc.end(), D[0].begin(), D[0].end());
        R.push_back(R.back()+D[0].size());
    }

    //Рассчитываем перестановку строк используя перестановку из
    //метода MaximumTransversal и метода приведения к блочному виду
    Pr.resize(N+1);
    for (size_t i = 1; i <= N; ++i) Pr[i] = PrMT[Pc[i]];

    //Рассчитываем транспонированную перестановку
    Prt.resize(N+1); Pct.resize(N+1);
    for (size_t i = 1; i <= N; ++i) Prt[Pr[i]] = Pct[Pc[i]] = i;

    //Заполняем структуру результата
    BBDS.swap(Pr, Pc, Prt, Pct, R, Is);

#ifdef TIME_LOG
#ifdef _OPENMP
    t1 = omp_get_wtime() - t;
    std::cout << "Рассчет перестановки " << t1 << std::endl;
    t = omp_get_wtime();
#endif
#endif
}


void MatrixBBDPreordering::MaximumTransversal(SparseMatrix &A, VectSizet &Pr)
{
    size_t N = A.H, N1 = N+1;

    Pr.resize(N1, 0);
    //for (size_t i = 1; i <= N; ++i) Pr[i] = i;
    //return;

    std::vector< std::vector<size_t> > H(N1);
    std::vector<size_t> visited(N1, 1), row_match(N1, 0), lookahead(N1, 0),
            colptrs(N1, 0), stack(N1, 0);

    size_t stack_col, col, row, t;

    for (size_t i = 1; i <= N; ++i)
    {
        for (size_t q = A.F[i]; q != SPARSE_END; q = A.N[q])
        {
            H[A.C[q]].push_back(i);
        }
    }

    size_t aug_no = 1, stack_last, stack_end;
    for (size_t i = 1; i <= N; ++i)
    {
        if (Pr[i] == 0) {
            stack[1] = i; stack_last = 1; stack_end = N1;
            colptrs[i] = 0;

            while(stack_last > 0) {
                stack_col = stack[stack_last];

                size_t j;
                for (j = lookahead[stack_col]; j < H[stack_col].size() &&
                     row_match[H[stack_col][j]] != 0; ++j){}

                lookahead[stack_col] = j + 1;

                if (j >= H[stack_col].size())
                {
                    size_t k;
                    for (k = colptrs[stack_col]; k < H[stack_col].size(); ++k)
                    {
                        t = visited[H[stack_col][k]];
                        if (t != aug_no && t != 0)
                            break;
                    }

                    colptrs[stack_col] = k + 1;

                    if (k >= H[stack_col].size())
                    {
                        --stack_last;
                        stack[--stack_end] = stack_col;
                        continue;
                    }
                    else
                    {
                        row = H[stack_col][k];
                        visited[row] = aug_no;
                        col = row_match[row];
                        stack[++stack_last] = col;
                        colptrs[col] = 0;
                    }
                }
                else
                {
                    row = H[stack_col][j];
                    visited[row] = aug_no;
                    while (row != 0)
                    {
                        col = stack[stack_last--];
                        t = Pr[col];
                        Pr[col] = row;
                        row_match[row] = col;
                        row = t;
                    }
                    ++aug_no;
                    break;
                }
            }

            if (Pr[i] == 0)
                for (size_t j = stack_end+1; j <= N; ++j)
                    visited[Pr[stack[j]]] = 0;
        }
    }
}

void MatrixBBDPreordering::SangiovannyVincentelli(
        VectVectSizet &H, VectVectSizet &D, VectSizet &Is,
        size_t Nmin, size_t Nmax)
{
    size_t Jsize, N = H.size()-1;
    VectSizet J;
    VectBool U(N+1,false), I(N+1,false);

    D.clear();
    D.push_back( VectSizet() ); //D_0
    D.push_back( VectSizet() ); //D_1

    Is.clear();
    Is.push_back(0);

    size_t m = 1; //Номер текущего блока

    //Шаг 1:
    //Найдем вершину с минимальным количеством связей
    size_t mu = 1, adjmin = ULONG_MAX, hs;
    for (size_t i = 1; i <= N; ++i)
    {
        if ((hs = H[i].size()) < adjmin)
        {
            adjmin = hs;
            mu = i;
        }
    }

    //И добавим её в D_1
    J = H[mu];
    Jsize = adjmin;
    D[m].push_back(mu);
    U[mu] = true;

    VectSizet optJ;
    size_t optJsize = ULONG_MAX, optdomainsize = 0, domainsize = 1;

    size_t unusedcount = N-1;
    while (unusedcount > 0)
    {
        if (Jsize == 0)
        {
            D.push_back( VectSizet() ); ++m;
            Is.push_back(Is.back());

            //Тот же Шаг 1, только ищем среди оставшихся
            mu = 1, adjmin = ULONG_MAX;
            for (size_t i = 1; i <= N; ++i)
            {
                if ((!U[i]) && ((hs = H[i].size()) < adjmin))
                {
                    adjmin = hs;
                    mu = i;
                }
            }

            J = H[mu];

            optJ.clear();
            optJsize = ULONG_MAX;
            optdomainsize = 0;
            domainsize = 0;
        }
        else
        {
            //Шаг 2:
            //Из разделителя выбираем элемент, который минимально увеличит
            //разделитель и добавляем его в текущий блок.

            VectSizet Jaddopt;
            size_t Jaddsizemin = ULONG_MAX;
            size_t muindex = 0;

            //Перебираем вершины разделителя и для каждой ищем вершины которые
            //добавятся в разделитель в случае добавления её в блок,
            //нас интересует минимальное количество новых в разделителе
            size_t i_end = J.size();
//#pragma omp parallel for
            for (size_t i = 0; i < i_end; ++i)
            {
                if (Jaddsizemin == 0) continue; //Костылик для OpenMP

                size_t cur = J[i];
                VectSizet Jadd;

                size_t adj, Jaddsize = 0;
                for (size_t k = 0; k < H[cur].size(); ++k)
                {
                    adj = H[cur][k];
                    if ((!U[adj]) &&
                        (std::find(J.begin(),J.end(), adj) == J.end()))
                    {
                        Jadd.push_back(adj);
                        if(++Jaddsize >= Jaddsizemin) break;
                    }
                }

                //Запоминаем новые вершины разделителя для "оптимальной" вершины
//#pragma omp critical
                if (Jaddsize < Jaddsizemin)
                {
                    mu = cur;
                    muindex = i;
                    Jaddsizemin = Jaddsize;
                    Jaddopt = Jadd;
                }

            }

            //Считаем новый разделитель
            J.erase(J.begin() + muindex);
            J.insert(J.end(), Jaddopt.begin(), Jaddopt.end());
        }

        Jsize = J.size(); //Запомним размер разделителя

        //Добавляем вершину в блок
        D[m].push_back(mu);
        U[mu] = true;
        --unusedcount;
        ++domainsize;

        //Шаг 3:
        //Определяем наимешьший разделитель при размере блока от Nmin до Nmax
        //Когда блок достигает максимального размера разделитель переносим в D_0,
        //а вершины что после последней вершины блока помечаем как еще не рассмотренные
        //(т.е. возвращаем из блока в исходный граф)
        if (domainsize >= Nmin)
        {
            //Пока <= Nmax запоминаем наименьший контур
            if (Jsize <= optJsize)
            {
                optJsize = Jsize;
                optJ = J;
                optdomainsize = domainsize;
            }

            if (domainsize == Nmax)
            {
                Is[m-1] += optJsize; //Запомним размер локального разделителя

                //Переносим вершины контура в интерфейс
                for (size_t j = 0; j < optJsize; ++j)
                {
                    size_t index = optJ[j];

                    //Сначала уберем связи данной вершины с графом
                    size_t k_end = H[index].size();
//#pragma omp parallel for
                    for (size_t k = 0; k < k_end; ++k)
                    {
                        size_t cur = H[index][k];
                        if (cur == index) continue;
                        VectSizetIterator it = std::find(H[cur].begin(),
                                                         H[cur].end(), index);
                        if (it != H[cur].end())
                            H[cur].erase(it);
                    }
                    H[index].clear();

                    //Пометим как вершину разделителя
                    I[index] = true;
                    if (!U[index])
                    {
                        U[index] = true;
                        --unusedcount;
                    }

                    //Добавим в общий разделитель
                    D[0].push_back(index);
                }

                //Вершины после последней вершины блока переносим в новый блок
                if (optdomainsize != domainsize)
                {
                    D.push_back( VectSizet() ); ++m;
                    Is.push_back(Is.back());

                    domainsize = 0;

                    //Переносим из текушего блока в новый блок,
                    //некоторые из них могут оказаться в разделителе, их не трогаем
                    for (VectSizetIterator it = D[m-1].begin() + optdomainsize,
                         it_end = D[m-1].end(); it != it_end; ++it)
                    {
                        if (!I[*it])
                        {
                            D[m].push_back(*it);
                            ++domainsize;
                        }
                    }

                    //Уберем лишние вершины из текущего блока
                    D[m-1].resize(optdomainsize);

                    //Правим текущий разделитель, убираем вершины которые уже в D_0.
                    for (VectSizetIterator it = J.begin(); it != J.end(); ++it)
                    {
                        if (I[*it] || U[*it])
                        {
                            J.erase(it); --it;
                        }
                    }

                    Jsize = J.size();

                    optJ.clear();
                    optJsize = ULONG_MAX;
                    optdomainsize = 0;
                }
                else
                {
                    Jsize = 0;
                }
            }
        }
    }
}

void MatrixBBDPreordering::SangiovannyVincentelliImproved(
        VectVectSizet &H, VectVectSizet &D, VectSizet &Is,
        size_t Nmin, size_t Nmax, size_t Dmax)
{
    size_t Jsize, N = H.size()-1;
    VectSizet J;
    VectBool U(N+1,false), I(N+1,false);

    D.clear();
    D.push_back( VectSizet() ); //D_0
    D.push_back( VectSizet() ); //D_1

    Is.clear();
    Is.push_back(0);

    size_t m = 1; //Номер текущего блока

    //Шаг 1:
    //Найдем вершину с минимальным количеством связей
    size_t mu = 1, adjmin = ULONG_MAX, hs;
    for (size_t i = 1; i <= N; ++i)
    {
        if ((hs = H[i].size()) < adjmin)
        {
            adjmin = hs;
            mu = i;
        }
    }

    //И добавим её в D_1
    J = H[mu];
    Jsize = adjmin;
    D[m].push_back(mu);
    U[mu] = true;

    VectSizet optJ;
    size_t optJsize = ULONG_MAX, optdomainsize = 0, domainsize = 1;

    size_t unusedcount = N-1;
    while (unusedcount > 0)
    {
        if (Jsize == 0)
        {
            D.push_back( VectSizet() ); ++m;
            Is.push_back(Is.back());

            //Тот же Шаг 1, только ищем среди оставшихся
            mu = 1, adjmin = ULONG_MAX;
            for (size_t i = 1; i <= N; ++i)
            {
                if ((!U[i]) && ((hs = H[i].size()) < adjmin))
                {
                    adjmin = hs;
                    mu = i;
                }
            }

            J = H[mu];

            optJ.clear();
            optJsize = ULONG_MAX;
            optdomainsize = 0;
            domainsize = 0;
        }
        else
        {
            //Шаг 2:
            //Из разделителя выбираем элемент, который минимально увеличит
            //разделитель и добавляем его в текущий блок.

            VectSizet Jaddopt;
            size_t Jaddsizemin = ULONG_MAX;
            size_t muindex = 0;

            //Перебираем вершины разделителя и для каждой ищем вершины которые
            //добавятся в разделитель в случае добавления её в блок,
            //нас интересует минимальное количество новых в разделителе
            size_t i_end = J.size();
#pragma omp parallel for
            for (size_t i = 0; i < i_end; ++i)
            {
                if (Jaddsizemin == 0) continue; //Костылик для OpenMP

                size_t cur = J[i];
                VectSizet Jadd;

                size_t adj, Jaddsize = 0;
                for (size_t k = 0; k < H[cur].size(); ++k)
                {
                    adj = H[cur][k];
                    if ((!U[adj]) &&
                        (std::find(J.begin(),J.end(), adj) == J.end()))
                    {
                        Jadd.push_back(adj);
                        if(++Jaddsize >= Jaddsizemin) break;
                    }
                }

                //Запоминаем новые вершины разделителя для "оптимальной" вершины
#pragma omp critical
                if (Jaddsize < Jaddsizemin)
                {
                    mu = cur;
                    muindex = i;
                    Jaddsizemin = Jaddsize;
                    Jaddopt = Jadd;
                }

            }

            //Считаем новый разделитель
            J.erase(J.begin() + muindex);
            J.insert(J.end(), Jaddopt.begin(), Jaddopt.end());
        }

        Jsize = J.size(); //Запомним размер разделителя

        //Добавляем вершину в блок
        D[m].push_back(mu);
        U[mu] = true;
        --unusedcount;
        ++domainsize;

        //Шаг 3:
        //Определяем наимешьший разделитель при размере блока от Nmin до Nmax
        //Когда блок достигает максимального размера разделитель переносим в D_0,
        //а вершины что после последней вершины блока помечаем как еще не рассмотренные
        //(т.е. возвращаем из блока в исходный граф)
        if (domainsize >= Nmin)
        {
            //Пока <= Nmax запоминаем наименьший контур
            if (Jsize <= optJsize)
            {
                optJsize = Jsize;
                optJ = J;
                optdomainsize = domainsize;
            }

            if (domainsize == Nmax)
            {
                Is[m-1] += optJsize; //Запомним размер локального разделителя

                //Переносим вершины контура в интерфейс
                for (size_t j = 0; j < optJsize; ++j)
                {
                    size_t index = optJ[j];

                    //Сначала уберем связи данной вершины с графом
                    size_t k_end = H[index].size();
#pragma omp parallel for
                    for (size_t k = 0; k < k_end; ++k)
                    {
                        size_t cur = H[index][k];
                        if (cur == index) continue;
                        VectSizetIterator it = std::find(H[cur].begin(),
                                                         H[cur].end(), index);
                        if (it != H[cur].end())
                            H[cur].erase(it);
                    }
                    H[index].clear();

                    //Пометим как вершину разделителя
                    I[index] = true;
                    if (!U[index])
                    {
                        U[index] = true;
                        --unusedcount;
                    }

                    //Добавим в общий разделитель
                    D[0].push_back(index);
                }

                //Вершины после последней вершины блока возвращаем в граф
                if (optdomainsize != domainsize)
                {
                    //Возвращаем в граф, помечая как не использованные
                    //некоторые из них могут оказаться в разделителе, их не трогаем
                    for (VectSizetIterator it = D[m].begin() + optdomainsize,
                         it_end = D[m].end(); it != it_end; ++it)
                    {
                        if (!I[*it])
                        {
                            U[*it] = false;
                            ++unusedcount;
                        }
                    }

                    //Уберем их из блока
                    D[m].resize(optdomainsize);
                }

                //Если размер разделителя достиг максимума, то неиспользованные
                //вершины переносим отдельный блок
                if (D[0].size() >= Dmax)
                {
                    D.push_back( VectSizet() ); ++m;
                    Is.push_back(Is.back());

                    for (size_t i = 1; i <= N; ++i)
                    {
                        if (!U[i]) D[m].push_back(i);
                    }

                    break; //Выходим, т.к. все вершины уже стали используемыми
                }

                Jsize = 0;
            }
        }
    }
}

#endif // !DECOMPOSITOR_H

