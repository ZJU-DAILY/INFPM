
// Created by xueqin on 2025/4/9.
//
#pragma once

#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <stack>
#include <random>
#include <algorithm>
#include "sfmt/SFMT.h"




using edge = Graph::edge;

class SCC {
public:

     sfmt_t sfmtSeed;
    // 构造函数：初始化 sfmtSeed



     void generate_scc_dags (
                    const TemporalGraph & tg, const Argument & arg,
                    unordered_map<int, unordered_map<int, vector<int>>>& FT_sccs, //T_sccs[t][scc_id]->访问第 scc_id 个 SCC 的所有节点
                    unordered_map<int, unordered_map<int, vector<pair<int, double>>>>& fin_all_scc_dags // all_scc_dags[t][r][scc_id] = list of (neighbor_scc_id, prob)
                    )
    {

        cout << "generate SCC DAG start. " << endl;
        int RSCC = arg._RSCC;
        int n = tg.size_n + 1;

        srand(time(NULL));
        sfmt_init_gen_rand(&sfmtSeed , time(NULL));

        int count = 0;
        for (int t : tg.T_valid) {
            //cout << "==== Time snapshot t = " << t << " ====" << std::endl;
            //cout << "count: " << count++ <<endl;

            vector<int> dfn(n, -1), low(n, -1), in_stack(n, 0), node_to_scc(n, -1);
            stack<int> stk;
            vector<edge> sampled_edges;

            int best_r = -1;
            int min_scc_cnt = INT_MAX;


            for (int r = 0; r < RSCC; ++r) {
                //std::cout << "-- RSCC sampling round r = " << r << " --" << std::endl;

                //// step 1: generate sampled graph
                sampled_edges.clear();
                for (const auto &e: tg.T_graph.at(t)) {
                    if (sfmt_genrand_real1(&sfmtSeed) < e.w) {
                        sampled_edges.emplace_back(e.st, e.ed, e.w);
                    }
                }


                // === Step 2: Tarjan's Algorithm to get SCCs ===
                unordered_map<int, vector<int>> sccs;
                fill(dfn.begin(), dfn.end(), -1);
                fill(low.begin(), low.end(), -1);
                fill(in_stack.begin(), in_stack.end(), 0);
                stack<int> empty;
                swap(stk, empty);
                //int idx = 0;
                //vector<int> node_to_scc(n, -1);
                fill(node_to_scc.begin(), node_to_scc.end(), -1);
                int idx = 0, scc_id = 0;

                std::function<void(int)> dfs = [&](int u) {
                    dfn[u] = low[u] = idx++;
                    stk.push(u);
                    in_stack[u] = 1;
                    for (const auto &e: sampled_edges) {
                        if (e.st != u) continue;
                        int v = e.ed;
                        if (dfn[v] == -1) {
                            dfs(v);
                            low[u] = std::min(low[u], low[v]);
                        } else if (in_stack[v]) {
                            low[u] = std::min(low[u], dfn[v]);
                        }
                    }
                    if (dfn[u] == low[u]) {
                        vector<int> scc;
                        while (true) {
                            int v = stk.top();
                            stk.pop();
                            in_stack[v] = 0;
                            scc.emplace_back(v);
                            node_to_scc[v] = scc_id;
                            if (v == u) break;
                        }
                        //sccs.emplace_back(std::move(scc));
                        sccs[scc_id] = std::move(scc);
                        scc_id++;
                    }
                };

                for (int i = 0; i < n; ++i) {
                    if (dfn[i] == -1) dfs(i);
                }


                if ((int) sccs.size() < min_scc_cnt) {
                    min_scc_cnt = sccs.size();
                    best_r = r;
                }
            }

            // === 重新执行 best_r 的采样 & SCC 提取（避免每轮保存）===
            sampled_edges.clear();
            for (const auto &e: tg.T_graph.at(t)) {
                if (sfmt_genrand_real1(&sfmtSeed) < e.w) {
                    sampled_edges.emplace_back(e.st, e.ed, e.w);
                }
            }

            // Final Tarjan for best_r
            unordered_map<int, vector<int>> best_sccs;
            fill(dfn.begin(), dfn.end(), -1);
            fill(low.begin(), low.end(), -1);
            fill(in_stack.begin(), in_stack.end(), 0);
            stack<int>().swap(stk);
            fill(node_to_scc.begin(), node_to_scc.end(), -1);
            int idx = 0, scc_id = 0;

            std::function<void(int)> dfs_final = [&](int u) {
                dfn[u] = low[u] = idx++;
                stk.push(u);
                in_stack[u] = 1;
                for (const auto &e: sampled_edges) {
                    if (e.st != u) continue;
                    int v = e.ed;
                    if (dfn[v] == -1) {
                        dfs_final(v);
                        low[u] = std::min(low[u], low[v]);
                    } else if (in_stack[v]) {
                        low[u] = std::min(low[u], dfn[v]);
                    }
                }
                if (dfn[u] == low[u]) {
                    vector<int> scc;
                    while (true) {
                        int v = stk.top();
                        stk.pop();
                        in_stack[v] = 0;
                        scc.emplace_back(v);
                        node_to_scc[v] = scc_id;
                        if (v == u) break;
                    }
                    best_sccs[scc_id] = std::move(scc);
                    scc_id++;
                }
            };

            for (int i = 0; i < n; ++i) {
                if (dfn[i] == -1) dfs_final(i);
            }

            // === Step 3: Build SCC-DAG and save ===
            unordered_map<int, int> scc_map;
            for (int v = 0; v < n; ++v) {
                if (node_to_scc[v] != -1) {
                    scc_map[v] = node_to_scc[v];
                }
            }

            FT_sccs[t] = std::move(best_sccs);

            fin_all_scc_dags[t] = build_scc_dag(tg.T_graph.at(t), scc_map);

        }


    }




    static unordered_map<int, std::vector<std::pair<int, double>>> build_scc_dag
                                (const vector<edge>& sampled_edges,
                                 const unordered_map<int, int>& node_to_scc) {

        std::unordered_map<long long, double> edge_prob;

        // Step 1: 合并跨 SCC 的边权（1 - ∏(1 - p)）
        for (const auto& e : sampled_edges) {
            int su = node_to_scc.at(e.st);
            int sv = node_to_scc.at(e.ed);
            if (su != sv) {
                long long key = ((long long)su << 32) | sv;
                edge_prob[key] = 1 - (1 - edge_prob[key]) * (1 - e.w);
            }
        }

        // Step 2: 构建带权的 SCC-DAG
        unordered_map<int, std::vector<std::pair<int, double>>> reverse_dag;
        for (const auto& kv : edge_prob) {
            int su = kv.first >> 32; //start
            int sv = kv.first & 0xffffffff; //end
            double prob = kv.second;

           // dag[su].emplace_back(sv, prob); // 添加边和边权
            reverse_dag[sv].emplace_back(su, prob);
        }

        return reverse_dag;

    }




};

