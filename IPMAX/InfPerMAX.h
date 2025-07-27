//
//
// Created by DBL-XQ on 2024/10/15.
//

#pragma once


class IPMAX
{
public:

    static void ReverseGreedy(TemporalGraph & tg, const Argument & arg)
    {
        cout << "ReverseGreedy begin..."<<endl;

        int k = int(tg.V * arg._k + 1);
        int theta = arg._theta;
        double alpha = arg._alpha;
        int size_n = tg.size_n;
        int R = arg._R;
        int T = tg.Tmax + 1;
        vector<int> & T_valid = tg.T_valid; // 将T_graph中包含的所有存在动作发生的时刻存入集合T_valid中
        cout << "valid time stamp num: " << T_valid.size() <<endl;


        double stime1 = omp_get_wtime();
        for (int i = 0; i < T_valid.size(); i++){
            int t_snap = T_valid[i];
            tg.buildhypergraph(R, t_snap);
        }
        double etime1 = omp_get_wtime();
        cout << "RR set generating time: " <<etime1 - stime1 <<endl;

        unordered_map<int, vector<int>> snapshot_seedsets; // 每个 snapshot 对应的 Si
        vector<bool> pruned(tg.Tmax+1, false);            // 是否被剪枝

        //// Step 1-1 选择每个 snapshot 下的 top-k 作为初始种子
        cout << "selected each snapshot graph top-k user by greedy" <<endl;
        for (int i : T_valid) {
            snapshot_seedsets[i] = tg.greedy_select_topk_nodes_by_rr(i, k, R);
        }

        //// Step 1-2 alpha 剪枝
        set<int> alpha_perst;
        double p = 0.5;
        for (int i : T_valid) {

            vector<int> can_seedset = snapshot_seedsets[i];
            vector<bool> rr_bitmap(R, false);
            for (int node : can_seedset) {
                for (int idx : tg.hyperG[i][node]) {
                    rr_bitmap[idx] = true;
                }
            }
            int coverTR = std::count(rr_bitmap.begin(), rr_bitmap.end(), true);
            double inf_t =(( double) coverTR / R) * (double) size_n;
            //cout << "inf_t: " << inf_t << " coverTR: " << coverTR <<endl;

            if (inf_t <= alpha * size_n * (1.0+p))  pruned[i] = true;
            else alpha_perst.insert(i);
        }

        cout << "[prune-alpha] total number of time stamps: " << alpha_perst.size()<<endl;



        //// Step 1-3 theta 剪枝（通过持久性）
        set<int> theta_perst;
        int init_perstence = 0;
        vector<pair<int, int>> persistent_inf_interval;
        if (alpha_perst.size() >= theta) {
            vector<int> seq(alpha_perst.begin(), alpha_perst.end());
            tg.findValidSequences(seq, theta, persistent_inf_interval);
            cout << "number of valid interval："<< persistent_inf_interval.size() <<endl;
            for (const auto& p : persistent_inf_interval) {
                int length = p.second - p.first +1;
                init_perstence += length;
                for (int t = p.first; t <=p.second; t++){
                    theta_perst.insert(t);
                }

            }
        }

        for (int t_snap : alpha_perst) {
            if (theta_perst.find(t_snap) == theta_perst.end()) {
                pruned[t_snap] = true;  // 在 alpha_perst 中但不在 theta_perst 中 → 被剪枝
            }
        }
        cout << "[prune-theta] total number of time stamps: " << theta_perst.size()<<endl;



        //// Step 1-4 生成候选集合 Cs （范围更小）
        unordered_set<int> candidate_set;         // 最终候选集合 Cs

        unordered_set<int> extra_nodes_set;
        for (int t : theta_perst) {
            for (int node : snapshot_seedsets[t])
                extra_nodes_set.insert(node);
        }
        vector<double> singnode_inf(size_n, 0.);
        for (unordered_set<int>::iterator it = extra_nodes_set.begin(); it != extra_nodes_set.end(); ++it) {
            int node = *it;
            for (int t_snap : theta_perst) {
                auto it_snap = tg.hyperG.find(t_snap);
                if (it_snap != tg.hyperG.end()) {
                    const auto& node_map = it_snap->second;
                    auto it_node = node_map.find(node);
                    if (it_node != node_map.end()) {
                        int coverTR = tg.hyperG[t_snap][node].size();
                        double inf_t = ((double) coverTR / R) * (double) size_n;
                        singnode_inf[node] += (double) inf_t;
                    }
                }

            }
        }

        vector<int> extra_nodes(extra_nodes_set.begin(), extra_nodes_set.end());
        sort(extra_nodes.begin(), extra_nodes.end(),
             [&singnode_inf](int a, int b) {
                 return singnode_inf[a] > singnode_inf[b];
             });
        for (int i = 0; i < extra_nodes.size(); ++i) {
            if (i < int(alpha * size_n) && singnode_inf[extra_nodes[i]]>0) candidate_set.insert(extra_nodes[i]);
        }



        cout << "candidate_set size: " << candidate_set.size() << endl ;
        cout << "selected seed size: " << k <<endl;

        /////////////////// step 2 从candidate_set中贪心删点得到最终的k个种子节点 //////////////////////

        vector<int> final_seedset(candidate_set.begin(), candidate_set.end());
        // 每个 snapshot 的 bitmap（RR 是否被覆盖）与覆盖次数
        set<int> loop_theta_perst = theta_perst;
        int current_persistence = init_perstence;
        vector<pair<int, int>> best_intervals;

        // 预处理每个节点覆盖的RR数量（覆盖数少的优先删）
        unordered_map<int, int> node_rr_coverage;
        for (int node : final_seedset) {
            int coverage = 0;
            for (int t : theta_perst) {
                coverage += tg.hyperG[t][node].size();
            }
            node_rr_coverage[node] = coverage;
        }

        // 贪心删点
        while ((int)final_seedset.size() > k) {
            int min_node = -1;
            int best_persistence = -1;
            int no_improve_cnt = 0;
            int no_improve_limit = 100;


            set<int> best_theta_perst;

            // 按node_rr_coverage升序考虑节点，优先删"贡献少"的
            vector<int> sorted_nodes(final_seedset.begin(), final_seedset.end());
            sort(sorted_nodes.begin(), sorted_nodes.end(), [&](int a, int b) {
                return node_rr_coverage[a] < node_rr_coverage[b];
            });

            //double stime1 = omp_get_wtime();
            for (int node : sorted_nodes) {
                if ( current_persistence - best_persistence == 0) break;

                vector<int> test_seedset;
                for (int x : final_seedset) {if (x != node) test_seedset.push_back(x);}

                set<int> loop_alpha_perst;
                for (int t : loop_theta_perst){
                    vector<bool> bitmap(R, false);
                    for (int node : test_seedset) {
                        for (int idx : tg.hyperG[t][node]) {
                            bitmap[idx] = true;
                        }
                    }
                    int cover = count(bitmap.begin(), bitmap.end(), true);
                    double inf_t = ((double)cover / R) * size_n;
                    if (inf_t >= alpha * size_n) {
                        loop_alpha_perst.insert(t);
                    }
                }
                int new_perst = 0;
                set<int> new_theta;
                vector<pair<int, int>> new_intervals;
                if ((int)loop_alpha_perst.size() >= theta) {
                    vector<int> seq(loop_alpha_perst.begin(), loop_alpha_perst.end());
                    tg.findValidSequences(seq, theta, new_intervals);
                    for (const auto& p : new_intervals) {
                        new_perst += p.second - p.first + 1;
                        for (int t = p.first; t <= p.second; ++t) {
                            new_theta.insert(t);
                        }
                    }
                }

                int drop = current_persistence - new_perst;
                if (min_node == -1 || drop < current_persistence - best_persistence) {
                    min_node = node;
                    best_persistence = new_perst;
                    best_theta_perst = new_theta;
                    best_intervals = new_intervals;
                    no_improve_cnt = 0;  // reset if improvement
                }else{
                    no_improve_cnt++;
                    if (no_improve_cnt >= no_improve_limit) {
                        //cout << "Early stop: no improvement in last " << no_improve_cnt << " attempts";
                        break;
                    }
                }
                //cout << "test node: " << node << " drop: " << drop << " best drop: " <<current_persistence - best_persistence<<endl;
            }


            if (min_node != -1) {
                final_seedset.erase(remove(final_seedset.begin(), final_seedset.end(), min_node), final_seedset.end());
                current_persistence = best_persistence;
                loop_theta_perst = best_theta_perst;
                persistent_inf_interval = best_intervals;



            } else {
                break;
            }
        }

        tg.seedset = final_seedset;



        cout << "=====================================================" <<endl;

        cout << "the max influence persistence: " << current_persistence <<endl;
        cout << "number of valid interval："<< persistent_inf_interval.size() <<endl;
        //cout << "candidate seedset size: " << candidate_set.size() << " candidate seeds: [ " ;
        //for (int can : candidate_set) {cout << can << ", ";}cout <<"]"<<endl;
        cout << "seedset size: " << tg.seedset.size() <<endl;
        cout << "seeds: [ " ;
        for (int val : tg.seedset) {
            cout << val << ", ";  // 输出每个值
        }cout <<"]"<<endl;
        //cout << "seeds: [ " ;for (int val : tg.seedset) {cout << val << ", ";}cout <<"]"<<endl;


    }


    static void ReverseGreedySCC(TemporalGraph & tg, const Argument & arg,
                                 const unordered_map<int, unordered_map<int,vector<int>>>& FT_sccs,
                                 const unordered_map<int, unordered_map<int, vector<pair<int, double>>>>& fin_all_scc_dags)
    {
        cout << "ReverseGreedySCC begin..."<< endl;

        int k = int(tg.V * arg._k + 1);
        int theta = arg._theta;
        double alpha = arg._alpha;
        int size_n = tg.size_n;
        int R = arg._R;
        int T = tg.Tmax + 1;
        vector<int> & T_valid = tg.T_valid; // 将T_graph中包含的所有存在动作发生的时刻存入集合T_valid中
        cout << "valid time stamp num: " << T_valid.size() <<endl;

        //// sample on SCC-DAG graph
        double stime1 = omp_get_wtime();
        for (int i = 0; i < T_valid.size(); i++){
            int t_snap = T_valid[i];
            tg.buildhypergraphSCC(R, t_snap, FT_sccs, fin_all_scc_dags);
        }
        double etime1 = omp_get_wtime();
        cout << "RR-SCC set generating time: " <<etime1 - stime1 <<endl;

        unordered_map<int, vector<int>> snapshot_seedsets; // 每个 snapshot 对应的 Si
        vector<bool> pruned(tg.Tmax+1, false);            // 是否被剪枝

        //// Step 1-1 选择每个 snapshot 下的 top-k 作为初始种子
        cout << "selected each snapshot graph top-k user by greedy" <<endl;
        for (int i : T_valid) {
            snapshot_seedsets[i] = tg.greedy_select_topk_nodes_by_rr(i, k, R);
        }

        //// Step 1-2 alpha 剪枝
        set<int> alpha_perst;
        double p = 0.5;
        for (int i : T_valid) {

            vector<int> can_seedset = snapshot_seedsets[i];
            vector<bool> rr_bitmap(R, false);
            for (int node : can_seedset) {
                for (int idx : tg.hyperG[i][node]) {
                    rr_bitmap[idx] = true;
                }
            }
            int coverTR = std::count(rr_bitmap.begin(), rr_bitmap.end(), true);
            double inf_t =(( double) coverTR / R) * (double) size_n;
            //cout << "inf_t: " << inf_t << " coverTR: " << coverTR <<endl;

            if (inf_t <= alpha * size_n * (1.0+p))  pruned[i] = true;
            else alpha_perst.insert(i);
        }

        cout << "[prune-alpha] total number of time stamps: " << alpha_perst.size()<<endl;


        //// Step 1-3 theta 剪枝（通过持久性）
        set<int> theta_perst;
        int init_perstence = 0;
        vector<pair<int, int>> persistent_inf_interval;
        if (alpha_perst.size() >= theta) {
            vector<int> seq(alpha_perst.begin(), alpha_perst.end());
            tg.findValidSequences(seq, theta, persistent_inf_interval);
            cout << "number of valid interval："<< persistent_inf_interval.size() <<endl;
            for (const auto& p : persistent_inf_interval) {
                int length = p.second - p.first +1;
                init_perstence += length;
                for (int t = p.first; t <=p.second; t++){
                    theta_perst.insert(t);
                }
                //cout << "Sequence (head: " << p.first << ", tail: " << p.second << ") "<<endl;
            }
        }

        for (int t_snap : alpha_perst) {
            if (theta_perst.find(t_snap) == theta_perst.end()) {
                pruned[t_snap] = true;  // 在 alpha_perst 中但不在 theta_perst 中 → 被剪枝
            }
        }
        cout << "[prune-theta] total number of time stamps: " << theta_perst.size()<<endl;



        //// Step 1-4 生成候选集合 Cs （范围更小）
        unordered_set<int> candidate_set;         // 最终候选集合 Cs

        unordered_set<int> extra_nodes_set;
        for (int t : theta_perst) {
            for (int node : snapshot_seedsets[t])
                extra_nodes_set.insert(node);
        }
        vector<double> singnode_inf(size_n, 0.);
        for (unordered_set<int>::iterator it = extra_nodes_set.begin(); it != extra_nodes_set.end(); ++it) {
            int node = *it;
            for (int t_snap : theta_perst) {
                auto it_snap = tg.hyperG.find(t_snap);
                if (it_snap != tg.hyperG.end()) {
                    const auto& node_map = it_snap->second;
                    auto it_node = node_map.find(node);
                    if (it_node != node_map.end()) {
                        int coverTR = tg.hyperG[t_snap][node].size();
                        double inf_t = ((double) coverTR / R) * (double) size_n;
                        singnode_inf[node] += (double) inf_t;
                    }
                }

            }
        }

        vector<int> extra_nodes(extra_nodes_set.begin(), extra_nodes_set.end());
        sort(extra_nodes.begin(), extra_nodes.end(),
             [&singnode_inf](int a, int b) {
                 return singnode_inf[a] > singnode_inf[b];
             });
        for (int i = 0; i < extra_nodes.size(); ++i) {
            if (i < int(alpha * size_n) && singnode_inf[extra_nodes[i]]>0) candidate_set.insert(extra_nodes[i]);
        }

        cout << "candidate_set size: " << candidate_set.size() << endl ;
        cout << "selected seed size: " << k <<endl;

        /////////////////// step 2 从candidate_set中贪心删点得到最终的k个种子节点 //////////////////////

        vector<int> final_seedset(candidate_set.begin(), candidate_set.end());
        // 每个 snapshot 的 bitmap（RR 是否被覆盖）与覆盖次数
        set<int> loop_theta_perst = theta_perst;
        int current_persistence = init_perstence;
        vector<pair<int, int>> best_intervals;

        // 预处理每个节点覆盖的RR数量（覆盖数少的优先删）
        unordered_map<int, int> node_rr_coverage;
        for (int node : final_seedset) {
            int coverage = 0;
            for (int t : theta_perst) {
                coverage += tg.hyperG[t][node].size();
            }
            node_rr_coverage[node] = coverage;
        }

        // 贪心删点
        while ((int)final_seedset.size() > k) {
            int min_node = -1;
            int best_persistence = -1;
            int no_improve_cnt = 0;
            int no_improve_limit = 100; //int(candidate_set.size()/2);


            set<int> best_theta_perst;

            // 按node_rr_coverage升序考虑节点，优先删"贡献少"的
            vector<int> sorted_nodes(final_seedset.begin(), final_seedset.end());
            sort(sorted_nodes.begin(), sorted_nodes.end(), [&](int a, int b) {
                return node_rr_coverage[a] < node_rr_coverage[b];
            });

            //double stime1 = omp_get_wtime();
            for (int node : sorted_nodes) {
                if ( current_persistence - best_persistence == 0) break;

                vector<int> test_seedset;
                for (int x : final_seedset) {if (x != node) test_seedset.push_back(x);}

                set<int> loop_alpha_perst;
                for (int t : loop_theta_perst){
                    vector<bool> bitmap(R, false);
                    for (int node : test_seedset) {
                        for (int idx : tg.hyperG[t][node]) {
                            bitmap[idx] = true;
                        }
                    }
                    int cover = count(bitmap.begin(), bitmap.end(), true);
                    double inf_t = ((double)cover / R) * size_n;
                    //double inf_t =((( double) cover / R)-sqrt(log(40.0) / R)) * (double) size_n;
                    if (inf_t >= alpha * size_n) {
                        loop_alpha_perst.insert(t);
                    }
                }
                int new_perst = 0;
                set<int> new_theta;
                vector<pair<int, int>> new_intervals;
                if ((int)loop_alpha_perst.size() >= theta) {
                    vector<int> seq(loop_alpha_perst.begin(), loop_alpha_perst.end());
                    tg.findValidSequences(seq, theta, new_intervals);
                    for (const auto& p : new_intervals) {
                        new_perst += p.second - p.first + 1;
                        for (int t = p.first; t <= p.second; ++t) {
                            new_theta.insert(t);
                        }
                    }
                }

                int drop = current_persistence - new_perst;
                if (min_node == -1 || drop < current_persistence - best_persistence) {
                    min_node = node;
                    best_persistence = new_perst;
                    best_theta_perst = new_theta;
                    best_intervals = new_intervals;
                    no_improve_cnt = 0;  // reset if improvement
                }else{
                    no_improve_cnt++;
                    if (no_improve_cnt >= no_improve_limit) {
                        //cout << "Early stop: no improvement in last " << no_improve_cnt << " attempts";
                        break;
                    }
                }
                //cout << "test node: " << node << " drop: " << drop << " best drop: " <<current_persistence - best_persistence<<endl;
            }


            if (min_node != -1) {
                final_seedset.erase(remove(final_seedset.begin(), final_seedset.end(), min_node), final_seedset.end());
                current_persistence = best_persistence;
                loop_theta_perst = best_theta_perst;
                persistent_inf_interval = best_intervals;


            } else {
                break;
            }
        }

        tg.seedset = final_seedset;



        cout << "=====================================================" <<endl;

        cout << "the max influence persistence: " << current_persistence <<endl;
        cout << "number of valid interval："<< persistent_inf_interval.size() <<endl;
        //cout << "candidate seedset size: " << candidate_set.size() << " candidate seeds: [ " ;
        //for (int can : candidate_set) {cout << can << ", ";}cout <<"]"<<endl;
        cout << "seedset size: " << tg.seedset.size() <<endl;
        cout << "seeds: [ " ;
        for (int val : tg.seedset) {
            cout << val << ", ";  // 输出每个值
        }cout <<"]"<<endl;
        //cout << "seeds: [ " ;for (int val : tg.seedset) {cout << val << ", ";}cout <<"]"<<endl;

    }



    static void LazyReplace_prn_opt(TemporalGraph & tg, const Argument & arg)
    {
        cout << "LazyReplace_prn begin..."<<endl;

        int k = int(tg.V * arg._k + 1);
        int theta = arg._theta;
        double alpha = arg._alpha;
        int size_n = tg.size_n;
        int R = arg._R;
        int T = tg.Tmax +1;
        vector<int> & T_valid = tg.T_valid; // 将T_graph中包含的所有存在动作发生的时刻存入集合T_valid中
        cout << "valid time stamp num: " << T_valid.size() <<endl;

        ////////////////////////////  1. 对每个t时刻snapshot图进行R次采样得到每个时刻下的采样集合   ////////////////////////
        double stime1 = omp_get_wtime();
        for (int i = 0; i < T_valid.size(); i++){
            int t_snap = T_valid[i];
            tg.buildhypergraph(R, t_snap);
        }
        double etime1 = omp_get_wtime();
        cout << "time1: " <<etime1 - stime1 <<endl;

        vector<int> can_seedset;
        vector<int> remain_nodes;
        //// 选择 inf-topk 和 freq-topk 两者中inf_persistence更大的作为初始种子集
        int inf_persistence = tg.find_best_initset(can_seedset, remain_nodes, k , R, alpha, theta);
        cout << "remain_set.size: " << remain_nodes.size() <<endl;

        int coverTR; double inf_t;
        vector<pair<int, int>> persistent_inf_interval;

        /////////////////////////////// start replace ////////////////////////////////////
        int count = 0;
        int best_inf_persistence = inf_persistence;  // 记录当前的最优影响力持久性
        set <int> best_perst;
        vector<int> best_seedset = can_seedset;  // 记录当前最优的种子集
        vector<pair<int, int>> best_persistent_inf_interval;
        unordered_map<int, int> best_seed_total_contrib;
        int v_weakest = -1;
        // 使用 remain_nodes 中的节点逐个替换 can_seedset 中的节点
        for (int x = 0; x < remain_nodes.size(); x++) {
            if (best_perst.size() == T_valid.size()) break; //已经达到所能影响到的最多有效时刻，退出替换
            int it = remain_nodes[x];

            if (x > 0){  //// prunning  x > 0
                // ========== [Step 2] 模拟替换 v_weakest 为 it ==========
                vector<int> temp_seed = best_seedset;
                replace(temp_seed.begin(), temp_seed.end(), v_weakest, it);
                set<int> curr_perst;
                for (int ti = 0; ti < T_valid.size(); ++ti) {
                    int t_snap = T_valid[ti];
                    vector<bool> rr_bitmap(R, false);
                    for (int node : temp_seed) {
                        for (int idx : tg.hyperG[t_snap][node]) {
                            rr_bitmap[idx] = true;
                        }
                    }
                    int coverTR = std::count(rr_bitmap.begin(), rr_bitmap.end(), true);
                    double inf_t = ((double)coverTR / R) * (double)size_n;
                    if (inf_t >= alpha * size_n) {
                        curr_perst.insert(t_snap);
                    }
                }

                // ========== [Step 3] 用 curr_perst 判断 interval 是否真的变长 ==========
                int new_inf_persistence = 0;
                if (curr_perst.size() >= theta) {
                    vector<int> seq(curr_perst.begin(), curr_perst.end());
                    vector<pair<int, int>> temp_intervals;
                    tg.findValidSequences(seq, theta, temp_intervals);
                    for (const auto& p : temp_intervals) {
                        new_inf_persistence += p.second - p.first + 1;
                    }
                }
                if (new_inf_persistence <= best_inf_persistence) {
                    //cout << "tight lower bound skip node: " << it << endl;
                    continue;  // 剪枝：替换最弱节点都无法提升，则此 it 可跳过
                }
            }

            //// not prunned + replace
            bool improved = false;
            int best_replaced_node = -1;  // 被替换的节点
            // 每次替换一个节点，替换完当前一轮之后，再选择 remain_nodes 中的下一个节点
            for (int i = can_seedset.size() - 1; i >= 0; --i) {
                double stime = omp_get_wtime();

                count ++;
                int replaced_node = can_seedset[i];  // 获取当前种子集中的节点
                can_seedset[i] = it;  // 将 remain_nodes 中的当前节点替换到种子集中的当前节点
                //cout << "replace node: " << can_seedset[i] << " replaced node: " << replaced_node <<endl;

                // 重新计算持久影响力时间段
                set<int> curr_perst;
                unordered_map<int, int> seed_total_contrib;
                for (int j = 0; j < T_valid.size(); j++){
                    int t_snap = T_valid[j];
                    //inf_t = tg.GetInfTgraph_Set_MC(R, seedset_cand, t_snap);
                    vector<bool> rr_bitmap(R, false);
                    for (int node : can_seedset) {
                        for (int idx : tg.hyperG[t_snap][node]) {
                            rr_bitmap[idx] = true;
                            seed_total_contrib[node]++;
                        }
                    }
                    coverTR = std::count(rr_bitmap.begin(), rr_bitmap.end(), true);
                    inf_t =(( double) coverTR / R) * (double) size_n;
                    if (inf_t >= alpha * size_n)
                        curr_perst.insert(t_snap);
                }
                vector<pair<int, int>> new_persistent_inf_interval;
                int new_inf_persistence = 0;  // 初始化新的持久影响力

                if (curr_perst.size() >= theta) {
                    vector<int> seq(curr_perst.begin(), curr_perst.end());
                    tg.findValidSequences(seq, theta, new_persistent_inf_interval);
                    //cout << "number of valid interval: " << new_persistent_inf_interval.size() << endl;
                    for (const auto& p : new_persistent_inf_interval) {
                        int length = p.second - p.first + 1;  // 计算长度
                        new_inf_persistence += length;     // 累加该节点所有序列的长度
                        //cout << "Sequence (head: " << p.first << ", tail: " << p.second << ")" << endl;
                    }

                    // 如果新的持久影响力大于当前最优的持久影响力，则更新最优种子集
                    if (new_inf_persistence > best_inf_persistence) {
                        best_perst = curr_perst;
                        best_inf_persistence = new_inf_persistence;
                        best_persistent_inf_interval = new_persistent_inf_interval;
                        best_seedset = can_seedset;
                        best_replaced_node = replaced_node;
                        best_seed_total_contrib = seed_total_contrib;
                        improved = true;
                    }
                }

                // 恢复种子集，准备进行下一个节点替换
                can_seedset[i] = replaced_node;  // 恢复之前的节点

                double etime = omp_get_wtime();
                //cout << "LR one round time: " << etime-stime<<endl;
            }

            // 如果有改进，更新种子集和持久影响力；否则退出循环
            if (improved) {
                // 更新最优种子集和持久影响力
                cout << "=========> Replace success, final replace node: " << it
                     << " replaced node: " << best_replaced_node << " count: " << count << endl;
                persistent_inf_interval = best_persistent_inf_interval;
                can_seedset = best_seedset;
                inf_persistence = best_inf_persistence;
                // ========== [Step 1] 预估当前最弱的种子节点 ==========
                int min_contrib = INT_MAX;
                for (const auto& pair : best_seed_total_contrib) {
                    int node = pair.first;
                    int score = pair.second;
                    if (score < min_contrib) {
                        v_weakest = node;
                        min_contrib = score;
                    }
                }

                cout << "=========> Replace success with the total sequence length: " << inf_persistence  << endl;
            }
        }


        tg.seedset = can_seedset;

        cout << "=====================================================" <<endl;
//        for (const auto& p : persistent_inf_interval) {
//            cout << "Sequence (head: " << p.first << ", tail: " << p.second << ") "<< "length: "<< p.second-p.first +1<<endl;
//        }
        cout << "the max influence persistence: " << inf_persistence <<endl;
        cout << "number of valid interval："<< persistent_inf_interval.size() <<endl;
        cout << "valid time stamp num: " << T_valid.size() <<endl;
        cout << "seedset size: " << tg.seedset.size() <<endl;
        cout << "seeds: [ " ;
        for (int val : tg.seedset) {
            cout << val << ", ";  // 输出每个值
        }cout <<"]"<<endl;

        //cout << "seeds: [ " ;for (int val : tg.seedset) {cout << val << ", ";}cout <<"]"<<endl;
        cout << "count: " << count<<endl;

    }


    static void LazyReplaceSCC(TemporalGraph & tg, const Argument & arg,
                               const unordered_map<int, unordered_map<int,vector<int>>>& FT_sccs,
                               const unordered_map<int, unordered_map<int, vector<pair<int, double>>>>& fin_all_scc_dags)
    {
        cout << "LazyReplaceSCC begin..."<< endl;

        int k = int(tg.V * arg._k + 1);
        int theta = arg._theta;
        double alpha = arg._alpha;
        int size_n = tg.size_n;
        int R = arg._R;
        int T = tg.Tmax +1;
        vector<int> & T_valid = tg.T_valid; // 将T_graph中包含的所有存在动作发生的时刻存入集合T_valid中
        cout << "valid time stamp num: " << T_valid.size() <<endl;

        ////////////////////////////  1. 对每个t时刻snapshot图进行R次采样得到每个时刻下的采样集合   ////////////////////////
        //// sample on SCC-DAG graph
        double stime1 = omp_get_wtime();
        for (int i = 0; i < T_valid.size(); i++){
            int t_snap = T_valid[i];
            //tg.buildhypergraph(R, t_snap);
            tg.buildhypergraphSCC(R, t_snap, FT_sccs, fin_all_scc_dags);
        }
        double etime1 = omp_get_wtime();
        cout << "RR-SCC set generating time: " <<etime1 - stime1 <<endl;

        vector<int> can_seedset;
        vector<int> remain_nodes;
        //// 选择 inf-topk 和 freq-topk 两者中inf_persistence更大的作为初始种子集
        int inf_persistence = tg.find_best_initset(can_seedset, remain_nodes, k , R, alpha, theta);
        cout << "remain_set.size: " << remain_nodes.size() <<endl;

        int coverTR; double inf_t;
        vector<pair<int, int>> persistent_inf_interval;

        /////////////////////////////// start replace ////////////////////////////////////
        int count = 0;
        int best_inf_persistence = inf_persistence;  // 记录当前的最优影响力持久性
        set <int> best_perst;
        vector<int> best_seedset = can_seedset;  // 记录当前最优的种子集
        vector<pair<int, int>> best_persistent_inf_interval;
        unordered_map<int, int> best_seed_total_contrib;
        int v_weakest = -1;
        double p=0.5;
        // 使用 remain_nodes 中的节点逐个替换 can_seedset 中的节点
        for (int x = 0; x < remain_nodes.size(); x++) {
            if (best_perst.size() == T_valid.size()) break; //已经达到所能影响到的最多有效时刻，退出替换
            int it = remain_nodes[x];

            if (x > 0){  //// prunning  x > 0
                // ========== [Step 2] 模拟替换 v_weakest 为 it ==========
                vector<int> temp_seed = best_seedset;
                replace(temp_seed.begin(), temp_seed.end(), v_weakest, it);
                set<int> curr_perst;
                for (int ti = 0; ti < T_valid.size(); ++ti) {
                    int t_snap = T_valid[ti];
                    vector<bool> rr_bitmap(R, false);
                    for (int node : temp_seed) {
                        for (int idx : tg.hyperG[t_snap][node]) {
                            rr_bitmap[idx] = true;
                        }
                    }
                    int coverTR = std::count(rr_bitmap.begin(), rr_bitmap.end(), true);
                    //double inf_t = ((double)coverTR / R) * (double)size_n;
                    double inf_t =((( double) coverTR / R)-sqrt(log(40.0) / R)) * (double) size_n;
                    if (inf_t >= alpha * size_n) {
                        curr_perst.insert(t_snap);
                    }
                }

                // ========== [Step 3] 用 curr_perst 判断 interval 是否真的变长 ==========
                int new_inf_persistence = 0;
                if (curr_perst.size() >= theta) {
                    vector<int> seq(curr_perst.begin(), curr_perst.end());
                    vector<pair<int, int>> temp_intervals;
                    tg.findValidSequences(seq, theta, temp_intervals);
                    for (const auto& p : temp_intervals) {
                        new_inf_persistence += p.second - p.first + 1;
                    }
                }
                if (new_inf_persistence <= best_inf_persistence) {
                    //cout << "tight lower bound skip node: " << it << endl;
                    continue;  // 剪枝：替换最弱节点都无法提升，则此 it 可跳过
                }
            }

            //// not prunned + replace
            bool improved = false;
            int best_replaced_node = -1;  // 被替换的节点
            // 每次替换一个节点，替换完当前一轮之后，再选择 remain_nodes 中的下一个节点
            for (int i = can_seedset.size() - 1; i >= 0; --i) {
                double stime = omp_get_wtime();

                count ++;
                int replaced_node = can_seedset[i];  // 获取当前种子集中的节点
                can_seedset[i] = it;  // 将 remain_nodes 中的当前节点替换到种子集中的当前节点
                //cout << "replace node: " << can_seedset[i] << " replaced node: " << replaced_node <<endl;

                // 重新计算持久影响力时间段
                set<int> curr_perst;
                unordered_map<int, int> seed_total_contrib;
                for (int j = 0; j < T_valid.size(); j++){
                    int t_snap = T_valid[j];
                    //inf_t = tg.GetInfTgraph_Set_MC(R, seedset_cand, t_snap);
                    vector<bool> rr_bitmap(R, false);
                    for (int node : can_seedset) {
                        for (int idx : tg.hyperG[t_snap][node]) {
                            rr_bitmap[idx] = true;
                            seed_total_contrib[node]++;
                        }
                    }
                    coverTR = std::count(rr_bitmap.begin(), rr_bitmap.end(), true);
                    //inf_t =(( double) coverTR / R) * (double) size_n;
                    inf_t =((( double) coverTR / R)-sqrt(log(40.0) / R)) * (double) size_n;
                    if (inf_t >=  alpha * size_n)
                        curr_perst.insert(t_snap);
                }
                vector<pair<int, int>> new_persistent_inf_interval;
                int new_inf_persistence = 0;  // 初始化新的持久影响力

                if (curr_perst.size() >= theta) {
                    vector<int> seq(curr_perst.begin(), curr_perst.end());
                    tg.findValidSequences(seq, theta, new_persistent_inf_interval);
                    //cout << "number of valid interval: " << new_persistent_inf_interval.size() << endl;
                    for (const auto& p : new_persistent_inf_interval) {
                        int length = p.second - p.first + 1;  // 计算长度
                        new_inf_persistence += length;     // 累加该节点所有序列的长度
                        //cout << "Sequence (head: " << p.first << ", tail: " << p.second << ")" << endl;
                    }

                    // 如果新的持久影响力大于当前最优的持久影响力，则更新最优种子集
                    if (new_inf_persistence > best_inf_persistence) {
                        best_perst = curr_perst;
                        best_inf_persistence = new_inf_persistence;
                        best_persistent_inf_interval = new_persistent_inf_interval;
                        best_seedset = can_seedset;
                        best_replaced_node = replaced_node;
                        best_seed_total_contrib = seed_total_contrib;
                        improved = true;
                    }
                }

                // 恢复种子集，准备进行下一个节点替换
                can_seedset[i] = replaced_node;  // 恢复之前的节点

                double etime = omp_get_wtime();

            }

            // 如果有改进，更新种子集和持久影响力；否则退出循环
            if (improved) {
                // 更新最优种子集和持久影响力
                cout << "=========> Replace success, final replace node: " << it
                     << " replaced node: " << best_replaced_node << " count: " << count << endl;
                persistent_inf_interval = best_persistent_inf_interval;
                can_seedset = best_seedset;
                inf_persistence = best_inf_persistence;
                // ========== [Step 1] 预估当前最弱的种子节点 ==========
                int min_contrib = INT_MAX;
                for (const auto& pair : best_seed_total_contrib) {
                    int node = pair.first;
                    int score = pair.second;
                    if (score < min_contrib) {
                        v_weakest = node;
                        min_contrib = score;
                    }
                }

                cout << "=========> Replace success with the total sequence length: " << inf_persistence  << endl;
            }
        }


        tg.seedset = can_seedset;

        cout << "=====================================================" <<endl;
//        for (const auto& p : persistent_inf_interval) {
//            cout << "Sequence (head: " << p.first << ", tail: " << p.second << ") "<< "length: "<< p.second-p.first +1<<endl;
//        }
        cout << "the max influence persistence: " << inf_persistence <<endl;
        cout << "number of valid interval："<< persistent_inf_interval.size() <<endl;
        cout << "valid time stamp num: " << T_valid.size() <<endl;
        cout << "seedset size: " << tg.seedset.size() <<endl;
        cout << "seeds: [ " ;
        for (int val : tg.seedset) {
            cout << val << ", ";  // 输出每个值
        }cout <<"]"<<endl;
        //cout << "seeds: [ " ;for (int val : tg.seedset) {cout << val << ", ";}cout <<"]"<<endl;
        cout << "count: " << count<<endl;

    }




//////////////////////////// window ////////////////////////////////////



    static void WinReverseGreedy(TemporalGraph & tg, const Argument & arg)
    {
        cout << "WinReverseGreedy begin..."<<endl;

        int k = int(tg.V * arg._k + 1);
        int theta = arg._theta;
        double alpha = arg._alpha;
        int size_n = tg.size_n;
        int R = arg._R;
        int W = arg._W;
        int T = tg.Tmax + 1;
        vector<int> & T_valid = tg.T_valid; // 将T_graph中包含的所有存在动作发生的时刻存入集合T_valid中
        cout << "valid time stamp num: " << T_valid.size() <<endl;
        //for (int t : T_valid) {std::cout << t << " ";}std::cout << std::endl;


        double stime1 = omp_get_wtime();
        for (int i = 0; i < T_valid.size(); i++){
            int t_snap = T_valid[i];
            tg.buildhypergraph(R, t_snap);
        }
        double etime1 = omp_get_wtime();
        cout << "RR set generating time: " <<etime1 - stime1 <<endl;

        unordered_map<int, vector<int>> snapshot_seedsets; // 每个 snapshot 对应的 Si
        vector<bool> pruned(tg.Tmax+1, false);            // 是否被剪枝

        //// Step 1-1 选择每个 snapshot 下的 top-k 作为初始种子
        cout << "selected each snapshot graph top-k user by greedy" <<endl;
        for (int i : T_valid) {
            snapshot_seedsets[i] = tg.greedy_select_topk_nodes_by_rr(i, k, R);
        }

        //// Step 1-2 alpha 剪枝
        set<int> alpha_perst;
        double p = 0.5;
        for (int i : T_valid) {

            vector<int> can_seedset = snapshot_seedsets[i];
            vector<bool> rr_bitmap(R, false);
            for (int node : can_seedset) {
                for (int idx : tg.hyperG[i][node]) {
                    rr_bitmap[idx] = true;
                }
            }
            int coverTR = std::count(rr_bitmap.begin(), rr_bitmap.end(), true);
            double inf_t =(( double) coverTR / R) * (double) size_n;
            //cout << "inf_t: " << inf_t << " coverTR: " << coverTR <<endl;

            if (inf_t <= alpha * size_n * (1.0+p))  pruned[i] = true;
            else alpha_perst.insert(i);
        }

        cout << "[prune-alpha] total number of time stamps: " << alpha_perst.size()<<endl;


        //// Step 1-3 theta 剪枝（通过持久性）
        set<int> theta_perst;
        int init_perstence = 0;
        vector<pair<int, int>> persistent_inf_interval;
        if (alpha_perst.size() >= theta) {
            vector<int> seq(alpha_perst.begin(), alpha_perst.end());
            // 查找所有连续子序列且满足尾 - 首 +1 >= timespan
            tg.findValidSequencesWin(seq, theta, W, persistent_inf_interval);
            // 计算每个有效序列的长度并累加
            cout << "number of valid interval："<< persistent_inf_interval.size() <<endl;
            for (const auto& p : persistent_inf_interval) {
                int length = p.second - p.first +1;
                init_perstence += length;
                for (int t = p.first; t <=p.second; t++){
                    theta_perst.insert(t);
                }

            }
        }

        for (int t_snap : alpha_perst) {
            if (theta_perst.find(t_snap) == theta_perst.end()) {
                pruned[t_snap] = true;  // 在 alpha_perst 中但不在 theta_perst 中 → 被剪枝
            }
        }

        if (theta_perst.size() > alpha_perst.size()) {
            theta_perst = alpha_perst;
        }
        cout << "[prune-theta] total number of time stamps: " << theta_perst.size()<<endl;



        //// Step 1-4 生成候选集合 Cs （范围更小）
        unordered_set<int> candidate_set;         // 最终候选集合 Cs

        unordered_set<int> extra_nodes_set;
        for (int t : theta_perst) {
            for (int node : snapshot_seedsets[t])
                extra_nodes_set.insert(node);
        }
        vector<double> singnode_inf(size_n, 0.);
        for (unordered_set<int>::iterator it = extra_nodes_set.begin(); it != extra_nodes_set.end(); ++it) {
            int node = *it;
            for (int t_snap : theta_perst) {
                auto it_snap = tg.hyperG.find(t_snap);
                if (it_snap != tg.hyperG.end()) {
                    const auto& node_map = it_snap->second;
                    auto it_node = node_map.find(node);
                    if (it_node != node_map.end()) {
                        int coverTR = tg.hyperG[t_snap][node].size();
                        double inf_t = ((double) coverTR / R) * (double) size_n;
                        singnode_inf[node] += (double) inf_t;
                    }
                }

            }
        }
        vector<int> extra_nodes(extra_nodes_set.begin(), extra_nodes_set.end());
        sort(extra_nodes.begin(), extra_nodes.end(),
             [&singnode_inf](int a, int b) {
                 return singnode_inf[a] > singnode_inf[b];
             });
        for (int i = 0; i < extra_nodes.size(); ++i) {
            if (i < int(alpha * size_n) && singnode_inf[extra_nodes[i]]>0) candidate_set.insert(extra_nodes[i]);
        }



        cout << "candidate_set size: " << candidate_set.size() << endl ;
        cout << "selected seed size: " << k <<endl;


        /////////////////// step 2 从candidate_set中贪心删点得到最终的k个种子节点 //////////////////////

        vector<int> final_seedset(candidate_set.begin(), candidate_set.end());
        // 每个 snapshot 的 bitmap（RR 是否被覆盖）与覆盖次数
        set<int> loop_theta_perst = theta_perst;
        int current_persistence = init_perstence;
        vector<pair<int, int>> best_intervals;

        // 预处理每个节点覆盖的RR数量（覆盖数少的优先删）
        unordered_map<int, int> node_rr_coverage;
        for (int node : final_seedset) {
            int coverage = 0;
            for (int t : theta_perst) {
                coverage += tg.hyperG[t][node].size();
            }
            node_rr_coverage[node] = coverage;
        }

        // 贪心删点
        while ((int)final_seedset.size() > k) {
            int min_node = -1;
            int best_persistence = -1;
            int no_improve_cnt = 0;
            int no_improve_limit = 100; //int(candidate_set.size()/2);

//          cout << "final_seedset size: " << final_seedset.size() << endl ;
//          cout << "final_seedset [ " ; for (int node : final_seedset) { cout << node << " ";}cout <<"]"<<endl;

            set<int> best_theta_perst;

            // 按node_rr_coverage升序考虑节点，优先删"贡献少"的
            vector<int> sorted_nodes(final_seedset.begin(), final_seedset.end());
            sort(sorted_nodes.begin(), sorted_nodes.end(), [&](int a, int b) {
                return node_rr_coverage[a] < node_rr_coverage[b];
            });

            //double stime1 = omp_get_wtime();
            for (int node : sorted_nodes) {
                if ( current_persistence - best_persistence == 0) break;

                vector<int> test_seedset;
                for (int x : final_seedset) {if (x != node) test_seedset.push_back(x);}

                set<int> loop_alpha_perst;
                for (int t : loop_theta_perst){
                    vector<bool> bitmap(R, false);
                    for (int node : test_seedset) {
                        for (int idx : tg.hyperG[t][node]) {
                            bitmap[idx] = true;
                        }
                    }
                    int cover = count(bitmap.begin(), bitmap.end(), true);
                    double inf_t = ((double)cover / R) * size_n;
                    if (inf_t >= alpha * size_n) {
                        loop_alpha_perst.insert(t);
                    }
                }
                int new_perst = 0;
                set<int> new_theta;
                vector<pair<int, int>> new_intervals;
                if ((int)loop_alpha_perst.size() >= theta) {
                    vector<int> seq(loop_alpha_perst.begin(), loop_alpha_perst.end());
                    tg.findValidSequencesWin(seq, theta, W, new_intervals);
                    for (const auto& p : new_intervals) {
                        new_perst += p.second - p.first + 1;
                        for (int t = p.first; t <= p.second; ++t) {
                            new_theta.insert(t);
                        }
                    }
                }

                int drop = current_persistence - new_perst;
                if (min_node == -1 || drop < current_persistence - best_persistence) {
                    min_node = node;
                    best_persistence = new_perst;
                    best_theta_perst = loop_alpha_perst; // 改为满足alpha条件的时刻数集合
                    best_intervals = new_intervals;
                    no_improve_cnt = 0;  // reset if improvement
                }else{
                    no_improve_cnt++;
                    if (no_improve_cnt >= no_improve_limit) {
                        //cout << "Early stop: no improvement in last " << no_improve_cnt << " attempts";
                        break;
                    }
                }
                //cout << "test node: " << node << " drop: " << drop << " best drop: " <<current_persistence - best_persistence<<endl;
            }
//            double etime1 = omp_get_wtime();
//            cout << "one round time: " << etime1 - stime1 << endl;

            if (min_node != -1) {
                final_seedset.erase(remove(final_seedset.begin(), final_seedset.end(), min_node), final_seedset.end());
                current_persistence = best_persistence;
                loop_theta_perst = best_theta_perst;
                persistent_inf_interval = best_intervals;

//                cout << "deleted [" << candidate_set.size()-final_seedset.size() << "] node: " << min_node
//                     << " current number of valid interval: " << persistent_inf_interval.size()
//                     << " current persistence: " << current_persistence <<endl;

            } else {
                break;
            }
        }

        tg.seedset = final_seedset;



        cout << "=====================================================" <<endl;
//        for (const auto& p : persistent_inf_interval) {
//            cout << "Sequence (head: " << p.first << ", tail: " << p.second << ") "<< "length: "<< p.second-p.first +1<<endl;
//        }
        cout << "the max influence persistence: " << current_persistence <<endl;
        cout << "number of valid interval："<< persistent_inf_interval.size() <<endl;
        //cout << "candidate seedset size: " << candidate_set.size() << " candidate seeds: [ " ;
        //for (int can : candidate_set) {cout << can << ", ";}cout <<"]"<<endl;
        cout << "seedset size: " << tg.seedset.size() <<endl;
        cout << "seeds: [ " ;
        for (int val : tg.seedset) {
            cout << val << ", ";  // 输出每个值
        }cout <<"]"<<endl;
        //cout << "seeds: [ " ;for (int val : tg.seedset) {cout << val << ", ";}cout <<"]"<<endl;


    }



    static void WinLazyReplace(TemporalGraph & tg, const Argument & arg)
    {
        cout << "LazyReplaceWin begin..."<<endl;

        int k = int(tg.V * arg._k + 1);
        int theta = arg._theta;
        double alpha = arg._alpha;
        int size_n = tg.size_n;
        int R = arg._R;
        int W = arg._W;
        int T = tg.Tmax + 1;
        vector<int> & T_valid = tg.T_valid; // 将T_graph中包含的所有存在动作发生的时刻存入集合T_valid中
        cout << "valid time stamp num: " << T_valid.size() <<endl;
        //for (int t : T_valid) {std::cout << t << " ";}std::cout << std::endl;

        ////////////////////////////  1. 对每个t时刻snapshot图进行R次采样得到每个时刻下的采样集合   ////////////////////////
        double stime1 = omp_get_wtime();
        for (int i = 0; i < T_valid.size(); i++){
            int t_snap = T_valid[i];
            tg.buildhypergraph(R, t_snap);
        }
        double etime1 = omp_get_wtime();
        cout << "time1: " <<etime1 - stime1 <<endl;

        vector<int> can_seedset;
        vector<int> remain_nodes;
        //// 选择 inf-topk 和 freq-topk 两者中inf_persistence更大的作为初始种子集
        int inf_persistence = tg.find_best_initset_win(can_seedset, remain_nodes, k , R, alpha, theta, W);
        cout << "remain_set.size: " << remain_nodes.size() <<endl;

        int coverTR; double inf_t;
        vector<pair<int, int>> persistent_inf_interval;


        /////////////////////////////// start replace ////////////////////////////////////
        int count = 0;
        int best_inf_persistence = inf_persistence;  // 记录当前的最优影响力持久性
        set <int> best_perst;
        vector<int> best_seedset = can_seedset;  // 记录当前最优的种子集
        vector<pair<int, int>> best_persistent_inf_interval = persistent_inf_interval;
        unordered_map<int, int> best_seed_total_contrib;
        int v_weakest = -1;
        // 使用 remain_nodes 中的节点逐个替换 can_seedset 中的节点
        for (int x = 0; x < remain_nodes.size(); x++) {
            if (best_perst.size() == T_valid.size()) break; //已经达到所能影响到的最多有效时刻，退出替换
            int it = remain_nodes[x];

            if (x > 0){  //// prunning  x > 0
                // ========== [Step 2] 模拟替换 v_weakest 为 it ==========
                vector<int> temp_seed = best_seedset;
                replace(temp_seed.begin(), temp_seed.end(), v_weakest, it);
                set<int> curr_perst;
                for (int ti = 0; ti < T_valid.size(); ++ti) {
                    int t_snap = T_valid[ti];
                    vector<bool> rr_bitmap(R, false);
                    for (int node : temp_seed) {
                        for (int idx : tg.hyperG[t_snap][node]) {
                            rr_bitmap[idx] = true;
                        }
                    }
                    int coverTR = std::count(rr_bitmap.begin(), rr_bitmap.end(), true);
                    double inf_t = ((double)coverTR / R) * (double)size_n;
                    if (inf_t >= alpha * size_n) {
                        curr_perst.insert(t_snap);
                    }
                }

                // ========== [Step 3] 用 curr_perst 判断 interval 是否真的变长 ==========
                int new_inf_persistence = 0;
                if (curr_perst.size() >= theta) {
                    vector<int> seq(curr_perst.begin(), curr_perst.end());
                    vector<pair<int, int>> temp_intervals;
                    tg.findValidSequencesWin(seq, theta, W, temp_intervals);
                    for (const auto& p : temp_intervals) {
                        new_inf_persistence += p.second - p.first + 1;
                    }
                }
                if (new_inf_persistence <= best_inf_persistence) {
                    //cout << "tight lower bound skip node: " << it << endl;
                    continue;  // 剪枝：替换最弱节点都无法提升，则此 it 可跳过
                }
            }

            //// not prunned + replace
            bool improved = false;
            int best_replaced_node = -1;  // 被替换的节点
            // 每次替换一个节点，替换完当前一轮之后，再选择 remain_nodes 中的下一个节点
            for (int i = can_seedset.size() - 1; i >= 0; --i) {
                // 替换当前的节点
                count ++;
                int replaced_node = can_seedset[i];  // 获取当前种子集中的节点
                can_seedset[i] = it;  // 将 remain_nodes 中的当前节点替换到种子集中的当前节点
                //cout << "replace node: " << can_seedset[i] << " replaced node: " << replaced_node <<endl;

                // 重新计算持久影响力时间段
                set<int> curr_perst;
                unordered_map<int, int> seed_total_contrib;
                for (int j = 0; j < T_valid.size(); j++){
                    int t_snap = T_valid[j];
                    //inf_t = tg.GetInfTgraph_Set_MC(R, seedset_cand, t_snap);
                    vector<bool> rr_bitmap(R, false);
                    for (int node : can_seedset) {
                        for (int idx : tg.hyperG[t_snap][node]) {
                            rr_bitmap[idx] = true;
                            seed_total_contrib[node]++;
                        }
                    }
                    coverTR = std::count(rr_bitmap.begin(), rr_bitmap.end(), true);
                    inf_t =(( double) coverTR / R) * (double) size_n;
                    if (inf_t >= alpha * size_n)
                        curr_perst.insert(t_snap);
                }


                // 如果 curr_perst 非空，计算当前种子集 can_seedset 的持久影响力时间段集合
                vector<pair<int, int>> new_persistent_inf_interval;
                int new_inf_persistence = 0;  // 初始化新的持久影响力

                if (curr_perst.size() >= theta) {
                    vector<int> seq(curr_perst.begin(), curr_perst.end());
                    tg.findValidSequencesWin(seq, theta, W,new_persistent_inf_interval);
                    //cout << "number of valid interval: " << new_persistent_inf_interval.size() << endl;
                    for (const auto& p : new_persistent_inf_interval) {
                        int length = p.second - p.first + 1;  // 计算长度
                        new_inf_persistence += length;     // 累加该节点所有序列的长度
                        //cout << "Sequence (head: " << p.first << ", tail: " << p.second << ")" << endl;
                    }


                    //cout << " with the total sequence length: " << new_inf_persistence << endl << endl;
                    // 如果新的持久影响力大于当前最优的持久影响力，则更新最优种子集
                    if (new_inf_persistence > best_inf_persistence) {
                        best_perst = curr_perst;
                        best_inf_persistence = new_inf_persistence;
                        best_persistent_inf_interval = new_persistent_inf_interval;
                        best_seedset = can_seedset;
                        best_replaced_node = replaced_node;
                        best_seed_total_contrib = seed_total_contrib;
                        improved = true;
                    }
                }

                // 恢复种子集，准备进行下一个节点替换
                can_seedset[i] = replaced_node;  // 恢复之前的节点
            }

            // 如果有改进，更新种子集和持久影响力；否则退出循环
            if (improved) {
                // 更新最优种子集和持久影响力
                cout << "=========> Replace success, final replace node: " << it
                     << " replaced node: " << best_replaced_node << " count: " << count << endl;
                persistent_inf_interval = best_persistent_inf_interval;
                can_seedset = best_seedset;
                inf_persistence = best_inf_persistence;
                // ========== [Step 1] 预估当前最弱的种子节点 ==========
                int min_contrib = INT_MAX;
                for (const auto& pair : best_seed_total_contrib) {
                    int node = pair.first;
                    int score = pair.second;
                    if (score < min_contrib) {
                        v_weakest = node;
                        min_contrib = score;
                    }
                }

                cout << "=========> Replace success with the total sequence length: " << inf_persistence  << endl;
            }
        }


        tg.seedset = can_seedset;

        cout << "=====================================================" <<endl;
//        for (const auto& p : persistent_inf_interval) {
//            cout << "Sequence (head: " << p.first << ", tail: " << p.second << ") "<< "length: "<< p.second-p.first +1<<endl;
//        }
        cout << "the max influence persistence: " << inf_persistence <<endl;
        cout << "number of valid interval："<< persistent_inf_interval.size() <<endl;
        cout << "valid time stamp num: " << T_valid.size() <<endl;
        cout << "seedset size: " << tg.seedset.size() <<endl;
        cout << "seeds: [ " ;
        for (int val : tg.seedset) {
            cout << val << ", ";  // 输出每个值
        }cout <<"]"<<endl;
        //cout << "seeds: [ " ;for (int val : tg.seedset) {cout << val << ", ";}cout <<"]"<<endl;
        cout << "count: " << count<<endl;

    }

};




