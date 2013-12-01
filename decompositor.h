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
    std::vector<size_t> G;
    size_t N;
    std::vector< Domain > D;

    std::vector<size_t> colorByBFS(size_t start, size_t color);
    void recolorGroup(Domain &Dt);
    void colorGroups();
    void moveNodeToInterface(size_t index, Domain &Dt);
    void moveNodeToInterface(size_t index);
    void decomposeMy();
    void decomposeDP(size_t M = 10);
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

    //decomposeMy();
    //decomposeDP();
    decomposeSangiovanniVincentelli();

    //Заполним структуру результата
    Pct.push_back(0);
    for (size_t i = 1; i < D.size(); ++i)
    {
        if (D[i].nodes.size() == 0) continue;
        Pct.insert(Pct.end(), D[i].nodes.begin(), D[i].nodes.end());
        R.push_back(D[i].nodes.size());
    }

    if (D[0].nodes.size() != 0)
    {
        Pct.insert(Pct.end(), D[0].nodes.begin(), D[0].nodes.end());
        R.push_back(D[0].nodes.size());
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
    G.clear();
    D.clear();
}

std::vector<size_t> MatrixDecompositor::colorByBFS(size_t start, size_t color)
{
    std::queue<size_t> nodes;
    std::vector<size_t> dnodes;
    size_t cur, t;

    nodes.push(start);
    G[start] = color;

    while (!nodes.empty())
    {
        cur = nodes.front();
        nodes.pop();

        dnodes.push_back(cur);

        for (size_t i = 0; i < H[cur].size(); ++i)
        {
            t = H[cur][i];
            if (G[t] == 0)
            {
                nodes.push(t);
                G[t] = color;
            }
        }
    }

    return dnodes;
}

void MatrixDecompositor::recolorGroup(Domain &Dt)
// Перекраска группы
{
    size_t N = Dt.nodes.size();
    size_t color = Dt.index;

    std::vector<size_t> nodesForRecol = Dt.nodes; //Запомним узлы, которые нужно перекрасить

    for (size_t k = 0; k < N; ++k)
        G[nodesForRecol[k]] = 0; //Сбросим индексы групп

    Dt.nodes = colorByBFS(nodesForRecol[0], color); //Раскрасим начиная с первого в группе

    color = D.size()-1;
    //Теперь если остались не раскрашенные, то раскрасим их
    size_t cur;
    for (size_t k = 1; k < N; ++k)
    {
        cur = nodesForRecol[k];
        if (G[cur] == 0)
        {
            ++color;
            D.push_back( Domain(color, colorByBFS(cur, color)) );
        }
    }
}

void MatrixDecompositor::colorGroups()
// Начальное разбиение на группы
{
    size_t color = 0; //Индекс текущей группы (группы будут начинаться с 1)
    G.assign(N, 0); //Сбросим индексы групп
    for (size_t k = 0; k < N; ++k)
    {
        if (G[k] == 0)
            //Если еще не занесен в группу
            //то делаем поиск в ширину и помечаем элементы текущим индексом
        {
            ++color;

            //std::vector<size_t> nodes = colorByBFS(k, color);
            D.push_back( Domain(color, colorByBFS(k, color)) );
        }
    }
}

void MatrixDecompositor::moveNodeToInterface(size_t index)
{
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
    G[index] = 0;
    D[0].nodes.push_back(index);
}

void MatrixDecompositor::moveNodeToInterface(size_t index, Domain& Dt)
{
    moveNodeToInterface(index);

    std::vector<size_t>::iterator it = std::find(Dt.nodes.begin(),
                                                 Dt.nodes.end(), index);
    if (it != Dt.nodes.end())
        Dt.nodes.erase(it);
}

void MatrixDecompositor::decomposeMy()
{
    std::vector< std::pair<size_t,size_t> > LS;  //Количество связей элементов
                                                 //в виде пары, чтобы можно было сортировать
    std::vector<bool> nodeFlags(N, true);
    D.clear();
    D.push_back( Domain(0, std::vector<size_t>()) );

    //Начальное определение групп
    colorGroups();

    bool nothingtoremove;

    while(true)
    {
        nothingtoremove = true;

        //Отсортируем группы по размеру по убыванию
        std::sort(D.begin()+1, D.end(), SortHelper::compareDomainSize);

        size_t M = D.size();

        double F = 0;
        for (size_t i = 1; i < M; ++i)
            F += pow(D[i].nodes.size(), 3.);
        F = pow(F, 1./3.);

        if (D[0].nodes.size() >= F) break;

        //Теперь последовательно пытаемся разделить самую большую группу
        for (size_t i = 1; i < M; ++i)
        {
            //Найдем количество связей
            LS.clear();
            for (size_t q = 0; q < D[i].nodes.size(); ++q)
            {
                size_t node = D[i].nodes[q];
                LS.push_back( std::pair<size_t,size_t>( node, H[node].size() ) );
            }

            //Отсортируем количество связей по убыванию
            std::sort(LS.begin(), LS.end(), SortHelper::comaprePairSecond);

            bool success = false;
            for (size_t j = 0; j < LS.size(); ++j)
            {
                //Теперь последовательно пробуем удалить элемент
                //с индексом LS[j].first
                //Он не должен оставлять элемент без связей
                //Если такого не найдется то выходим
                //
                size_t index = LS[j].first;
                if (!nodeFlags[index]) continue;

                //Проверим на то, что он не оставит элементов без связей
                bool goodNode = true;
                for (size_t k = 0; k < H[index].size(); ++k)
                    if (H[ H[index][k] ].size() <= 1)
                    {
                        goodNode = false;
                        break;
                    }

                if (!goodNode)
                {
                    nodeFlags[index] = false;
                    continue;
                }

                //Убираем элемент
                moveNodeToInterface(index, D[i]);
                //И перекрашиваем группу, вдруг разделилась
                recolorGroup(D[i]);
                success = true;
                break;
            }

            if (success)
            {
                nothingtoremove = false;
                break;
            }
        }

        if (nothingtoremove) break;
    }
}

void MatrixDecompositor::decomposeDP(size_t M)
{
    D.clear();
    D.push_back( Domain(0, std::vector<size_t>()) );

    G.assign(N, 0);

    //Начальное определение групп
    //colorGroups();

    srand(time(NULL));

    std::vector<size_t> V;

    M = std::max(2.,(log10(N*1.) - 1.) * 10.);

    std::vector< std::pair<size_t,size_t> > LS;

    for (size_t i = 0; i < H.size(); ++i)
        LS.push_back( std::pair<size_t,size_t>( i, H[i].size() ) );

    std::sort(LS.begin(), LS.end(), SortHelper::comaprePairSecondInv);

    for (size_t i = 0; i < M; ++i)
    {
        //size_t q = rand()%(N/M) + i*(N/M);
        size_t q = LS[i].first;
        V.push_back(q);
        G[q] = i+1;
        D.push_back( Domain(i+1, std::vector<size_t>(1, q)) );
    }

    size_t nNodes = M;

    //while (nNodes < N)
    //{
        for (size_t vi = 0; vi < V.size(); ++vi)
        {
            size_t i = V[vi];
            size_t color = G[i];

            if (color == M+1) continue; //Вершина из интерфейса

            for (size_t j = 0; j < H[i].size(); ++j) //Красим соседей
            {
                size_t cur = H[i][j];

                if (G[cur] == color) continue; //Если цвет уже такой же, идем дальше
                if (G[cur] == M+1) continue; //Вершина из интерфейса, идем дальше

                if (G[cur] == 0) //Если не раскрашена, красим в текущий цвет
                {
                    G[cur] = color; //Красим
                    D[color].nodes.push_back(cur); //Добавляем в группу
                    ++nNodes;
                    V.push_back(cur);
                }
                else if (G[cur] != M+1) //Если помечена другим цветом, то переносим в интефейс
                {
                    D[0].nodes.push_back(cur);

                    for (size_t k = 0; k < D[G[cur]].nodes.size(); ++k)
                        if (D[G[cur]].nodes[k] == cur)
                        {
                            D[G[cur]].nodes.erase(D[G[cur]].nodes.begin()+k);
                            break;
                        }

                    G[cur] = M+1;
                }
            }
        }
    //}
    D.push_back( Domain(M+1, std::vector<size_t>()) );
    for (size_t i = 0; i < N; ++i)
        if (G[i] == 0)
            D[M+1].nodes.push_back(i);
}

void MatrixDecompositor::decomposeSangiovanniVincentelli()
{
    std::vector<size_t> cn;
    std::vector<size_t> is;
    std::vector<size_t> asn, asp;
    std::vector<bool> isInterface(N, false);

    G.assign(N,N);

    D.clear();
    D.push_back( Domain(0, std::vector<size_t>()) );
    D.push_back( Domain(1, std::vector<size_t>()) );

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
    D.back().nodes.push_back(in);
    G[in] = 1;

    for (size_t i = 1; i < N; ++i)
    {
        asp = asn; //Запомним AS(i-1), вдруг это bottleneck? Тогда будем переносить в интерфейс.

        //Что делать если CN(i-1) == 0?
        //Делаем так же как и на первом шаге для тех элементов, что не в группах
        //и не в интерфейсе
        if (cn[i-1] == 0)
        {
            in = 0, asmin = ULONG_MAX;
            for (size_t i = 0; i < N; ++i)
            {
                if ((G[i] == N) && (!isInterface[i]) && (H[i].size() < asmin))
                {
                    asmin = H[i].size();
                    in = i;
                }
            }
            asn = H[in];
            D.push_back( Domain(1, std::vector<size_t>()) );
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
                    if ((std::find(asp.begin(),asp.end(), adj) == asp.end()) &&
                        (G[adj] == N))
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
        if ((D.back().nodes.size() > 0.75*Nmax) && (cn[i-2] >= cn[i-1]) && (cn[i-1] < cn[i]))
        {
            for (size_t j = 0; j < asp.size(); ++j)
            {
                size_t rem = asp[j];
                if (isInterface[rem]) continue;
                moveNodeToInterface(rem);
                isInterface[rem] = true;
            }

            D.push_back( Domain(1, std::vector<size_t>()) );
        }
        else
        {
            if (!isInterface[in])
            {
                D.back().nodes.push_back(in);
                G[in] = 1;
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
