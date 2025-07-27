//
// Created by DBL-XQ on 2025/1/7.
//

#pragma once


#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cassert>
#include <ctime>
#include <cmath>

#include <string>
#include <memory>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>
#include <list>
#include <limits>
#include <algorithm>
#include <map>
#include <unordered_map>
#include <set>
#include <stack>
#include <queue>
#include <omp.h>
#include <unordered_set>
#include <unordered_map>
#include "sfmt/SFMT.h"


using namespace std;


class Graph
{
public:
    int V, E; // the original node size and edge size
    int size_n = 0; // the max node_id+1 in graph



    int Tmin = INT32_MAX;
    int Tmax;

    int rangeTime = 86400; // 1 day =86400


    set<int> nodeset;
    //unordered_map<int, int> isTarget;

    //vector<vector<pair<int, int>>> T_graph;

    struct edge {
        int st;      // 起始节点
        int ed;      // 结束节点
        double w;   // 邻接节点的权重
        edge() {}
        edge(int a, int b, double c) : st(a), ed(b), w(c) {}
    };

    unordered_map<int, vector<edge>> T_graph;
    unordered_map<int, unordered_map<int, int>> inDeg;




    sfmt_t sfmtSeed;

    string attribute_file;
    string graph_file;
    int r;

    void readNM()
    {
        std::ifstream infile(attribute_file);
        if (!infile.is_open())
        {
            std::cout << "The file \"" + attribute_file + "\" can NOT be opened\n";
            exit(255);
            return;
        }

        infile >> V >> E;
        cout << "V: " << V << "\tE: " <<E<< endl;
        infile.close();

        //cout << "Read attribute end." << endl;

    }


    void readGraph(){

        std::ifstream infile(graph_file);
        if (!infile.is_open())
        {
            std::cout << "The file \"" + graph_file + "\" can NOT be opened\n";
            exit(255);
            return;
        }

        int timespan = 0;
        for (int i = 0; i < E; ++i)
        {
            int u, v, t;
            infile >> u >> v >> t;

            if (u == v) continue;
            if (u >= V || v>= V) continue;

            int t_group = t/rangeTime; // compute time interval

            if (t_group > timespan) timespan = t_group;
            if (t_group < Tmin) Tmin = t_group;

//            if (t_group >= T_graph.size()) {
//                T_graph.resize(t_group + 1);
//            }
            T_graph[t_group].emplace_back(u, v, 0.);

            nodeset.insert(u);
            nodeset.insert(v);
            if (u >= size_n) size_n = u;
            if (v >= size_n) size_n = v;
            inDeg[t_group][v]++;  // 更新归类时间段的入度
        }
        infile.close();

        size_n++;
        cout << "size_n: "<< size_n<<endl;


        Tmax = timespan;
        cout << "Tmin: "<< Tmin << ", Tmax: "<< Tmax<< endl;



        //计算并赋予权重 1/d_in,注意此处每个时刻的权重都不同
        //#pragma omp parallel for num_threads(4) reduction(+:time1, time2)
        //#pragma omp parallel for reduction(+:time1, time2)

        for (auto it = T_graph.begin(); it != T_graph.end(); ++it) {
            int tt = it->first;
            const auto& edges = it->second;
            if (T_graph[tt].size() == 0)
                continue;

            // 遍历 T_graph[t] 中的每个边 (u, v)
            for (const auto& edge : T_graph[tt]) {
                int u = edge.st;   // 起点
                int v = edge.ed;  // 终点

                // 判断终点
                double weight;

                //double startTime2 = omp_get_wtime();
                if (inDeg[tt][v] == 0) {weight = 0.;}
                else{
                    weight = 1.0 / inDeg[tt][v];
                    // 更新边的权重
                    for (auto& e : T_graph[tt]) {
                        if (e.ed == v) {
                            e.w = weight;  // 更新权重
                        }
                    }
                }
            }
        }


        //cout << "Assign weights end." <<endl;
        cout << "Read graph end." << endl;


    }


    Graph(string attribute_file, string graph_file): attribute_file(attribute_file), graph_file(graph_file)
    {
        srand(time(NULL));
        sfmt_init_gen_rand(&sfmtSeed, rand());

        readNM();

        readGraph();

    }

};


















