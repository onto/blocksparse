#ifndef DECOMPOSITOR_H
#define DECOMPOSITOR_H

#include <cstdlib>
#include <ctime>
#include <vector>
#include <algorithm>
#include <queue>
#include <cstring>

#include "sparsematrix.h"
#include "blockmatrix.h"

struct DecompositorResult
{
    DecompositorResult() {}

    std::vector<size_t> P; //Перестановка
    std::vector<size_t> R; //Размерности блоков
};

struct Domain
{
    Domain(size_t index_, std::vector<size_t> nodes_):
        index(index_),nodes(nodes_){}

    size_t index;
    std::vector<size_t> nodes;
};

class SortHelper
{
public:
    static bool comaprePairSecond(const std::pair<size_t,size_t>& p1, const std::pair<size_t,size_t>& p2)
    {
        return (p1.second > p2.second);
    }

    static bool compareDomainSize(const Domain& d1, const Domain& d2)
    {
        return (d1.nodes.size() > d2.nodes.size());
    }
};

class MatrixDecompositor
{
public:
    MatrixDecompositor(const char *file);
    ~MatrixDecompositor() {}

    DecompositorResult getResult();

private:
    std::vector< std::vector<size_t> > H;
    std::vector<size_t> G;
    size_t N;
    std::vector< Domain > D;

    DecompositorResult res;

    std::vector<size_t> colorByBFS(size_t start, size_t color);
    void recolorGroup(size_t index);
    void colorGroups();
    void decomposeMy();
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
                H[i].push_back(c-1);
            } else {
                break;
            }
        }
    }
    in.close();

    //std::cout << std::endl;
    //printH();

    decomposeMy();

    //Освободим память
    H.clear();
    G.clear();
    D.clear();
}

DecompositorResult MatrixDecompositor::getResult()
{
    return res;
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

void MatrixDecompositor::recolorGroup(size_t index)
// Перекраска группы
{
    size_t N = D[index].nodes.size();
    size_t color = D[index].index;

    for (size_t k = 0; k < N; ++k)
        G[D[index].nodes[k]] = 0; //Сбросим индексы групп

    std::vector<size_t> group1 = colorByBFS(D[index].nodes[0], color); //Раскрасим начиная с первого в группе

    color = D.size();

    //Теперь если остались не раскрашенные, то раскрасим их
    size_t cur;
    for (size_t k = 1; k < N; ++k)
    {
        cur = D[index].nodes[k];
        if (G[cur] == 0)
        {
            ++color;

            std::vector<size_t> nodes = colorByBFS(cur, color);
            D.push_back( Domain(color, nodes) );
        }
    }

    D[index].nodes = group1;
}

void MatrixDecompositor::colorGroups()
// Начальное разбиение на группы
{
    size_t color = 0; //Индекс текущей группы (группы будут начинаться с 1)
    size_t N = G.size();
    G.assign(N, 0); //Сбросим индексы групп
    for (size_t k = 0; k < N; ++k)
    {
        if (G[k] == 0)
            //Если еще не занесен в группу
            //то делаем поиск в ширину и помечаем элементы текущим индексом
        {
            ++color;

            std::vector<size_t> nodes = colorByBFS(k, color);
            D.push_back( Domain(color, nodes) );
        }
    }
}

void MatrixDecompositor::decomposeMy()
{
    G.assign(N, N);
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

        if (D[0].nodes.size() >= D[1].nodes.size()) break;

        //Теперь последовательно пытаемся разделить самую большую группу
        size_t M = D.size();
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
            for (size_t j = 0; j < N; ++j)
            {
                //Теперь последовательно пробуем удалить элемент
                //с индексом LS[j].first
                //Он не должен оставлять элемент без связей
                //Если такого не найдется то выходим
                //
                size_t index = LS[j].first;
                if (!nodeFlags[index]) continue;

//                //Проверка на то, что элемент не останется без связей
//                if (H[index].size() <= 1)
//                {
//                    nodeFlags[index] = false;
//                    continue;
//                }

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
                for (size_t k = 0; k < H[index].size(); ++k)
                {
                    size_t cur = H[index][k];
                    if (cur == index) continue;
                    for (size_t l = 0; l < H[cur].size(); ++l)
                    {
                        if (H[cur][l] == index)
                        {
                            H[cur].erase(H[cur].begin() + l);
                            break;
                        }
                    }
                }
                H[index].clear();
                G[index] = 0;
                for (size_t k = 1; k < D[i].nodes.size(); ++k)
                    if (D[i].nodes[k] == index) 
                    {
                        D[i].nodes.erase(D[i].nodes.begin()+k);
                        break;
                    }
                D[0].nodes.push_back(index);

                success = true;
                break;
            }

            if (success)
            {
                recolorGroup(i);
                nothingtoremove = false;
                break;
            }
        }

        if (nothingtoremove) break;
    }

    //Заполним структуру результата
    for (size_t i = 1; i < D.size(); ++i)
    {
        res.P.insert(res.P.end(), D[i].nodes.begin(), D[i].nodes.end());
        res.R.push_back(D[i].nodes.size());
    }

    res.P.insert(res.P.end(), D[0].nodes.begin(), D[0].nodes.end());
    res.R.push_back(D[0].nodes.size());
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
