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
#include "udecompositor.h"

class SortHelper
{
public:
    static bool comaprePairSecond(const std::pair<size_t,size_t>& p1,
                                  const std::pair<size_t,size_t>& p2)
    {
        return (p1.second > p2.second);
    }

    static bool comaprePairSecondInv(const std::pair<size_t,size_t>& p1,
                                     const std::pair<size_t,size_t>& p2)
    {
        return (p1.second < p2.second);
    }

    static bool compareDomainSize(const Domain& d1, const Domain& d2)
    {
        return (d1.nodes.size() > d2.nodes.size());
    }
};

class MatrixDecompositor : public SimpleMatrixDecompositor
{
public:
    MatrixDecompositor(const char *file);
    ~MatrixDecompositor() {}

private:
    std::vector< std::vector<size_t> > H;
    size_t N;
    std::vector< std::vector<size_t> > D;
    std::vector<bool> U, I;

    void moveNodeToInterface(size_t index);
    void decomposeSangiovanniVincentelli();
    void printH();
};

MatrixDecompositor::MatrixDecompositor(const char *file)
{
    std::ifstream in(file);

    in >> N >> N;

    H.assign(N, std::vector<size_t>());

    double x;
    size_t c;

    //Сначала считаем только связи
    for (size_t r = 1, i = 0; r <= N; ++r, ++i)
    {
        while (true)
        {
            in >> c;

            if (c != SPARSE_END){
                in >> x;
                if (r != c) H[i].push_back(c-1);
            } else {
                break;
            }
        }
    }
    in.close();

    decomposeSangiovanniVincentelli();

    //Заполним структуру результата
    Pct.push_back(0);
    for (size_t i = 1; i < D.size(); ++i)
    {
        if (D[i].size() == 0) continue;
        Pct.insert(Pct.end(), D[i].begin(), D[i].end());
        R.push_back(D[i].size());
    }

    if (D[0].size() != 0)
    {
        Pct.insert(Pct.end(), D[0].begin(), D[0].end());
        R.push_back(D[0].size());
    }

    Pr.resize(N+1);
    for (size_t i = 1; i <= N; ++i)
    {
        Pct[i] += 1;
        Pr[Pct[i]] = i;
    }

    std::ofstream out("permut.txt");
    for (size_t i = 0; i < Pr.size(); ++i)
        out << Pr[i] << " ";
    out.close();

    for (size_t i = 0; i < R.size(); ++i)
        std::cout << R[i] << " ";
    std::cout << std::endl;


    //Освободим память
    H.clear();
    D.clear();
    U.clear();
    I.clear();
}

void MatrixDecompositor::moveNodeToInterface(size_t index)
{
#pragma omp parallel for
    for (size_t k = 0; k < H[index].size(); ++k)
    {
        size_t cur = H[index][k];
        if (cur == index) continue;
        std::vector<size_t>::iterator it = std::find(H[cur].begin(),
                                                     H[cur].end(), index);
        if (it != H[cur].end())
            H[cur].erase(it);
    }
    H[index].clear();
    I[index] = true; U[index] = true;
    D[0].push_back(index);
}

void MatrixDecompositor::decomposeSangiovanniVincentelli()
{
    std::vector<size_t> cn;
    std::vector<size_t> is;
    std::vector<size_t> asn, asp;

    U.assign(N,false);
    I.assign(N,false);

    D.clear();
    D.push_back( std::vector<size_t>() );
    D.push_back( std::vector<size_t>() );

    //Вот это определим сами
    size_t Nmax = N / 10.;

    //Step 1
    //Найдем initial iterating node
    //с минимальным количеством связей
    size_t in = 0, asmin = ULONG_MAX;
    for (size_t i = 0; i < N; ++i)
    {
        if (H[i].size() < asmin)
        {
            asmin = H[i].size();
            in = i;
        }
    }

    //Step 2-4
    asn = H[in];
    cn.push_back(asmin);
    is.push_back(in);
    D.back().push_back(in);
    U[in] = true;

    for (size_t i = 1; i < N; ++i)
    {
        asp = asn; //Запомним AS(i-1), вдруг это bottleneck? Тогда будем переносить в интерфейс.

        //Что делать если CN(i-1) == 0?
        //Делаем так же как и на первом шаге для тех элементов, что не в группах
        if (cn[i-1] == 0)
        {
            in = 0, asmin = ULONG_MAX;
            for (size_t i = 0; i < N; ++i)
            {
                if ((!U[i]) && (H[i].size() < asmin))
                {
                    asmin = H[i].size();
                    in = i;
                }
            }
            asn = H[in];
            D.push_back( std::vector<size_t>() );
        }
        else
        {
            //Step 6
            //Выбор next iterating node
            //используя greedly strategy

            //Мы должны выбрать IS(i) с минимальным CN(i)
            std::vector<size_t> asadd, asaddopt;
            size_t addmin = ULONG_MAX;
            size_t inindx = 0;

#pragma omp parallel for private (asadd)
            for (size_t j = 0; j < asp.size(); ++j)
            {
                size_t cur = asp[j];
                asadd.clear();

                //Переберем соседей, если они не в IS(0) U ... U IS(i-1), те
                //соседей которые еще не входят ни в одну группу
                //и не в AS(i-1), то добавим в asadd
                for (size_t k = 0; k < H[cur].size(); ++k)
                {
                    size_t adj = H[cur][k];
                    if ((!U[adj]) &&
                        (std::find(asp.begin(),asp.end(), adj) == asp.end()))
                    {
                        asadd.push_back(adj);
                    }
                }

#pragma omp critical
                if (asadd.size() < addmin)
                {
                    in = cur;
                    inindx = j;
                    addmin = asadd.size();
                    asaddopt = asadd;
                }
            }
            asn.erase(asn.begin()+inindx);
            asn.insert(asn.end(),asaddopt.begin(),asaddopt.end());
        }

        //Step 7-8
        cn.push_back(asn.size());
        is.push_back(in);

        //Вроде как определение локального минимума
        //И тут перенос в интерфейс всех узлов из AS(i-1)
        if ((D.back().size() > 0.75*Nmax) && (cn[i-2] >= cn[i-1]) && (cn[i-1] < cn[i]))
        {
            for (size_t j = 0; j < asp.size(); ++j)
            {
                size_t rem = asp[j];
                if (I[rem]) continue;
                moveNodeToInterface(rem);
            }

            D.push_back( std::vector<size_t>() );
        }
        else
        {
            if (!I[in])
            {
                D.back().push_back(in);
                U[in] = true;
            }
        }
    }
}


void MatrixDecompositor::printH()
{
    for (size_t i = 0; i < H.size(); ++i)
    {
        for (size_t j = 0; j < H[i].size(); ++j)
            std::cout << H[i][j] << "\t";
        std::cout << std::endl;
    }
}

#endif // !DECOMPOSITOR_H
