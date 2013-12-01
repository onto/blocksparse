#ifndef UDECOMPOSITOR_H
#define UDECOMPOSITOR_H

#include <cstdlib>
#include <ctime>
#include <vector>
#include <algorithm>
#include <queue>
#include <cstring>
#include <climits>

#include "sparsematrix.h"
#include "blockmatrix.h"

struct Domain
{
    Domain(size_t index_, std::vector<size_t> nodes_):
        index(index_),nodes(nodes_){}

    size_t index;
    std::vector<size_t> nodes;
};

class SimpleMatrixDecompositor
{
public:
    SimpleMatrixDecompositor() {}

    std::vector<size_t> Pr, Prt, Pc, Pct, R;
};

class MatrixUDecompositor : public SimpleMatrixDecompositor
{
public:
    MatrixUDecompositor(const char *file);
    ~MatrixUDecompositor() {}

private:
    std::vector< std::vector<size_t> > Hr;
    std::vector< std::vector<size_t> > Hc;
    std::vector<size_t> Gr;
    std::vector<size_t> Gc;
    std::vector< Domain > Dr;
    std::vector< Domain > Dc;
    std::vector<bool> isInterface_r, isInterface_c;

    size_t N;

    void moveNodeToInterface_r(size_t index);
    void moveNodeToInterface_c(size_t index);
    void decomposeSangiovanniVincentelli();
    //void printH();
};

MatrixUDecompositor::MatrixUDecompositor(const char *file):
    SimpleMatrixDecompositor()
{
    std::ifstream in(file);

    in >> N >> N;

    Hr.assign(N, std::vector<size_t>());
    Hc.assign(N, std::vector<size_t>());

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
                if (r != c)
                {
                    Hr[i].push_back(c-1);
                    Hc[c-1].push_back(r-1);
                }
            } else {
                break;
            }
        }
    }
    in.close();

    //decomposeMy();
    //decomposeDP();
    decomposeSangiovanniVincentelli();

    //Заполним структуру результата
    Prt.push_back(0); Pct.push_back(0);
    for (size_t i = 1; i < Dr.size(); ++i)
    {
        if (Dr[i].nodes.size() == 0) continue;
        Prt.insert(Prt.end(), Dr[i].nodes.begin(), Dr[i].nodes.end());
        Pct.insert(Pct.end(), Dc[i].nodes.begin(), Dc[i].nodes.end());
        R.push_back(Dr[i].nodes.size());
    }

    if (Dr[0].nodes.size() != 0)
    {
        Prt.insert(Prt.end(), Dr[0].nodes.begin(), Dr[0].nodes.end());
        Pct.insert(Pct.end(), Dc[0].nodes.begin(), Dc[0].nodes.end());
        R.push_back(Dr[0].nodes.size());
    }

    Pr.resize(N+1); Pc.resize(N+1);
    for (size_t i = 1; i <= N; ++i)
    {
        Prt[i] += 1;
        Pr[Prt[i]] = i;
        Pct[i] += 1;
        Pc[Pct[i]] = i;
    }

    std::ofstream out_r("permut_r.txt");
    std::ofstream out_c("permut_c.txt");
    for (size_t i = 0; i <= N; ++i)
    {
        out_r << Pr[i] << " ";
        out_c << Pc[i] << " ";
    }
    out_r.close(); out_c.close();

    for (size_t i = 0; i < R.size(); ++i)
        std::cout << R[i] << " ";
    std::cout << std::endl;


    //Освободим память
    Hr.clear();Hr.clear();
    Dr.clear();Dc.clear();
    Gr.clear();Gc.clear();
}

void MatrixUDecompositor::moveNodeToInterface_r(size_t index)
{
    if (isInterface_r[index]) return;
    for (size_t k = 0; k < Hr[index].size(); ++k)
    {
        size_t cur = Hr[index][k];
        if (cur == index) continue;
        std::vector<size_t>::iterator it = std::find(Hr[cur].begin(),
                                                     Hr[cur].end(), index);
        if (it != Hr[cur].end())
            Hr[cur].erase(it);
    }
    Hr[index].clear();
    Gr[index] = 0;
    Dr[0].nodes.push_back(index);
    isInterface_r[index] = true;
}

void MatrixUDecompositor::moveNodeToInterface_c(size_t index)
{
    if (isInterface_c[index]) return;
    for (size_t k = 0; k < Hc[index].size(); ++k)
    {
        size_t cur = Hc[index][k];
        if (cur == index) continue;
        std::vector<size_t>::iterator it = std::find(Hc[cur].begin(),
                                                     Hc[cur].end(), index);
        if (it != Hc[cur].end())
            Hc[cur].erase(it);
    }
    Hc[index].clear();
    Gc[index] = 0;
    Dc[0].nodes.push_back(index);
    isInterface_c[index] = true;
}

void MatrixUDecompositor::decomposeSangiovanniVincentelli()
{
    std::vector<size_t> cn_r, cn_c;
    std::vector<size_t> is_r, is_c;
    std::vector<size_t> asn_r, asn_c, asp_r, asp_c;

    isInterface_r.assign(N, false);
    isInterface_c.assign(N, false);

    Gr.assign(N,N);
    Gc.assign(N,N);

    Dr.clear(); Dc.clear();
    Dr.push_back( Domain(0, std::vector<size_t>()) );
    Dc.push_back( Domain(0, std::vector<size_t>()) );
    Dr.push_back( Domain(1, std::vector<size_t>()) );
    Dc.push_back( Domain(1, std::vector<size_t>()) );

    //Вот это определим сами
    size_t Nmax = N / 10.;

    //Найдем initial iterating node для строк
    //с минимальным количеством связей
    size_t in_r = 0, in_c = 0, asmin = ULONG_MAX;
    for (size_t i = 0; i < N; ++i)
    {
        if (Hr[i].size() < asmin)
        {
            asmin = Hr[i].size();
            in_r = i;
        }
    }

    asn_r = Hr[in_r];
    cn_r.push_back(asmin);
    is_r.push_back(in_r);
    Dr.back().nodes.push_back(in_r);
    Gr[in_c] = 1;

    //Теперь для столбцов
    asmin = ULONG_MAX;
    for (size_t i = 0; i < N; ++i)
    {
        if (Hc[i].size() < asmin)
        {
            asmin = Hc[i].size();
            in_c = i;
        }
    }

    asn_c = Hc[in_c];
    cn_c.push_back(asmin);
    is_c.push_back(in_c);
    Dc.back().nodes.push_back(in_c);
    Gc[in_c] = 1;

    for (size_t i = 1; i < N; ++i)
    {
        asp_r = asn_r; //Запомним AS(i-1), вдруг это bottleneck?
        asp_c = asn_c; //Тогда будем переносить в интерфейс.

        //Что делать если CN(i-1) == 0?
        //Делаем так же как и на первом шаге для тех элементов, что не в группах
        //и не в интерфейсе
        if ((cn_r[i-1] == 0) || (cn_c[i-1] == 0))
        {
            if (cn_r[i-1] == 0)
            {
                in_r = 0, asmin = ULONG_MAX;
                for (size_t i = 0; i < N; ++i)
                {
                    if ((Gr[i] == N) && (!isInterface_r[i]) && (Hr[i].size() < asmin))
                    {
                        asmin = Hr[i].size();
                        in_r = i;
                    }
                }
                asn_r = Hr[in_r];
            }
            else
            {
                for (size_t j = 0; j < asp_r.size(); ++j)
                {
                    moveNodeToInterface_r(asp_r[j]);
                }
            }

            if (cn_c[i-1] == 0)
            {
                in_c = 0, asmin = ULONG_MAX;
                for (size_t i = 0; i < N; ++i)
                {
                    if ((Gc[i] == N) && (!isInterface_c[i]) && (Hc[i].size() < asmin))
                    {
                        asmin = Hc[i].size();
                        in_c = i;
                    }
                }
                asn_c = Hc[in_c];
            }
            else
            {
                for (size_t j = 0; j < asp_c.size(); ++j)
                {
                    moveNodeToInterface_r(asp_c[j]);
                }
            }

            Dr.push_back( Domain(1, std::vector<size_t>()) );
            Dc.push_back( Domain(1, std::vector<size_t>()) );
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

            for (size_t j = 0; j < asp_r.size(); ++j)
            {
                size_t cur = asp_r[j];

                if (isInterface_r[cur]) continue;

                asadd.clear();

                //Переберем соседей, если они не в IS(0) U ... U IS(i-1), те
                //соседей которые еще не входят ни в одну группу
                //и не в AS(i-1), то добавим в asadd
                for (size_t k = 0; k < Hr[cur].size(); ++k)
                {
                    size_t adj = Hr[cur][k];
                    if ((std::find(asp_r.begin(),asp_r.end(), adj) == asp_r.end()) &&
                        (Gr[adj] == N))
                    {
                        asadd.push_back(adj);
                    }

                }

                if (asadd.size() < addmin)
                {
                    in_r = cur;
                    inindx = j;
                    addmin = asadd.size();
                    asaddopt = asadd;
                }
            }
            asn_r.erase(asn_r.begin()+inindx);
            asn_r.insert(asn_r.end(),asaddopt.begin(),asaddopt.end());


            addmin = ULONG_MAX; inindx = 0;
            asaddopt.clear();
            for (size_t j = 0; j < asp_c.size(); ++j)
            {
                size_t cur = asp_c[j];

                if (isInterface_c[cur]) continue;

                asadd.clear();

                //Переберем соседей, если они не в IS(0) U ... U IS(i-1), те
                //соседей которые еще не входят ни в одну группу
                //и не в AS(i-1), то добавим в asadd
                for (size_t k = 0; k < Hc[cur].size(); ++k)
                {
                    size_t adj = Hc[cur][k];
                    if ((std::find(asp_c.begin(),asp_c.end(), adj) == asp_c.end()) &&
                        (Gc[adj] == N))
                    {
                        asadd.push_back(adj);
                    }

                }

                if (asadd.size() < addmin)
                {
                    in_c = cur;
                    inindx = j;
                    addmin = asadd.size();
                    asaddopt = asadd;
                }
            }
            asn_c.erase(asn_c.begin()+inindx);
            asn_c.insert(asn_c.end(),asaddopt.begin(),asaddopt.end());
        }

        //Step 7-8
        cn_r.push_back(asn_r.size());
        is_r.push_back(in_r);
        cn_c.push_back(asn_c.size());
        is_c.push_back(in_c);


        //Вроде как определение локального минимума
        //И тут перенос в интерфейс всех узлов из AS(i-1)
        if ((Dr.back().nodes.size() > 0.75*Nmax) &&
                (((cn_r[i-2] >= cn_r[i-1]) && (cn_r[i-1] < cn_r[i])) ||
                 ((cn_c[i-2] >= cn_c[i-1]) && (cn_c[i-1] < cn_c[i]))) )
        {
            for (size_t j = 0; j < asp_r.size(); ++j)
            {
                moveNodeToInterface_r(asp_r[j]);
            }
            for (size_t j = 0; j < asp_c.size(); ++j)
            {
                moveNodeToInterface_c(asp_c[j]);
            }

            Dr.push_back( Domain(1, std::vector<size_t>()) );
            Dc.push_back( Domain(1, std::vector<size_t>()) );
        }
        else
        {
            if (!isInterface_r[in_r])
            {
                Dr.back().nodes.push_back(in_r);
                Gr[in_r] = 1;
            }
            if (!isInterface_c[in_c])
            {
                Dc.back().nodes.push_back(in_c);
                Gc[in_c] = 1;
            }
        }
    }
}


//void MatrixUDecompositor::printH()
//{
//    for (size_t i = 0; i < H.size(); ++i)
//    {
//        for (size_t j = 0; j < H[i].size(); ++j)
//            std::cout << H[i][j] << "\t";
//        std::cout << std::endl;
//    }
//}




#endif // UDECOMPOSITOR_H
