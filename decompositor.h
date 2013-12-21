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

class MatrixDecompositor : public SimpleMatrixDecompositor
{
public:
    MatrixDecompositor(const char *file);
    ~MatrixDecompositor() {}

private:
    std::vector< std::vector<size_t> > H;
    size_t N, unusedcount;
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
    if (I[index]) return;

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
    I[index] = true;
    if (!U[index])
    {
        U[index] = true;
        --unusedcount;
    }
    D[0].push_back(index);
}

void MatrixDecompositor::decomposeSangiovanniVincentelli()
{
    std::vector<size_t> cn;
    std::vector<size_t> asn;

    U.assign(N,false);
    I.assign(N,false);

    D.clear();
    D.push_back( std::vector<size_t>() );
    D.push_back( std::vector<size_t>() );

    //Вот это определим сами
    size_t Nmax = std::max(2, int(N) / 10);
    size_t Nmin = 0.75*Nmax;

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
    D.back().push_back(in);
    U[in] = true;

    std::vector<size_t> mincontour;
    size_t minconsize = ULONG_MAX, dsize = 0;

    unusedcount = N-1;
    size_t i = 0;
    while (unusedcount > 0)//for (size_t i = 1; i < N; ++i)
    {
        ++i;

        //Что делать если CN(i-1) == 0?
        //Делаем так же как и на первом шаге для тех элементов, что не в группах
        if (cn[i-1] == 0)
        {
            in = 0, asmin = ULONG_MAX;
            for (size_t j = 0; j < N; ++j)
            {
                if ((!U[j]) && (H[j].size() < asmin))
                {
                    asmin = H[j].size();
                    in = j;
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
            for (size_t j = 0; j < asn.size(); ++j)
            {
                size_t cur = asn[j];
                asadd.clear();

                //Переберем соседей, если они не в IS(0) U ... U IS(i-1), те
                //соседей которые еще не входят ни в одну группу
                //и не в AS(i-1), то добавим в asadd
                size_t adj, addsize = 0;
                for (size_t k = 0; k < H[cur].size(); ++k)
                {
                    adj = H[cur][k];
                    if ((!U[adj]) &&
                        (std::find(asn.begin(),asn.end(), adj) == asn.end()))
                    {
                        asadd.push_back(adj);
                        if(++addsize >= addmin) break;
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

        if (!I[in])
        {
            D.back().push_back(in);
            U[in] = true;
            --unusedcount;
        }

        if (D.back().size() >= Nmin)
        {
            if (cn.back() <= minconsize)
            {
                minconsize = cn.back();
                mincontour = asn;
                dsize = D.back().size();
            }

            if (D.back().size() == Nmax)
            {
                for (size_t j = 0; j < mincontour.size(); ++j)
                    moveNodeToInterface(mincontour[j]);

                if (dsize == D.back().size())
                {
                    cn.back() = 0;
                }
                else
                {
                    std::vector<size_t> dadd;
                    dadd.assign(D.back().begin()+dsize, D.back().end());
                    D.back().resize(dsize);

                    D.push_back( std::vector<size_t>() );

                    for (size_t j = 0; j < dadd.size(); ++j)
                    {
                        if (!I[dadd[j]])
                            D.back().push_back(dadd[j]);

                        std::vector<size_t>::iterator it = std::find(asn.begin(), asn.end(), dadd[j]);
                        if (it != asn.end())
                            asn.erase(it);
                    }

                    if (D.back().size() == 0)
                    {
                        D.pop_back();
                        cn.back() = 0;
                    }
                }

                minconsize = ULONG_MAX;
                mincontour.clear();
                dsize = 0;
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
