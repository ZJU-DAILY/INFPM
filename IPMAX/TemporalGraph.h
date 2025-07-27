
//
// Created by DBL-XQ on 2025/3/9.
//


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


#include "graph.h"

class TemporalGraph: public Graph
{
private:
    unordered_map<int, vector<bool>>visit;
    unordered_map<int, vector<int>> visit_mark;



public:
    int validTimeNum;
    vector<int> T_valid;
    int avg_size_n = 0; // the avg number of nodes in graph + 1
    unordered_map<int, set<int>> targetsOfT_graph; // target node set of each snapshot graph
    vector<int> seedset; // the final seedset

    //// TR set
    //unordered_map<int, vector<vector<int>>> hyperG; //  hyperG[t][u] = {rr_id1, rr_id2, ...}
    unordered_map<int, vector<vector<int>>> hyperGT; // hyperGT[t][rr_id] = {u1, u2, ...}

    unordered_map<int, unordered_map<int, vector<int>>> hyperG;
    //unordered_map<int, unordered_map<int, vector<int>>> hyperGT;

    sfmt_t sfmtSeed;

    TemporalGraph(string attribute_file, string graph_file) : Graph(attribute_file, graph_file)
    {

        srand(time(NULL));
        sfmt_init_gen_rand(&sfmtSeed , time(NULL));
        //sfmt_init_gen_rand(&sfmtSeed, 95082); //By using a determined seed number, we could debug without randomness.

        validTimeNum = 0;
        int total_active_node_count = 0;
        for (const auto& it : T_graph) {
            int t = it.first;
            const auto& edges = it.second;

            if (!edges.empty()) {
                validTimeNum++;
                T_valid.push_back(t);

                visit[t] = vector<bool>(size_n, false);
                visit_mark[t] = vector<int>(size_n);
                hyperGT[t] = vector<vector<int>>();

                unordered_set<int> active_nodes;  // 当前 snapshot 的活跃节点
                for (const auto& e : edges) {
                    active_nodes.insert(e.st);
                    active_nodes.insert(e.ed);
                }
                total_active_node_count += active_nodes.size();
            }
        }

        avg_size_n = validTimeNum > 0 ? int((double)total_active_node_count / validTimeNum) : 0;

        cout << "total valid time num: " << validTimeNum << endl;
        cout << "avg active node num: " << avg_size_n << endl;
        std::cout << "Initial Temporal graph end. " << endl;


    }



    void buildhypergraph(unsigned int R, const int timeStamp) // build TR sets of graph at timeStamp
    {
        if (R > INT_MAX){
            cout << "Error:R too large" << endl;
            exit(1);
        }
        //ASSERT(hyperGT.size() == 0);
        unsigned int size = hyperGT[timeStamp].size();
        int uStart;
        while (size < R){
            uStart = sfmt_genrand_uint32(&sfmtSeed) % (size_n);
            BuildHypergraphNode(uStart, size++, timeStamp);
        }

    }

    void BuildHypergraphNode(const int uStart, unsigned int hyperId, const int timeStamp)
    {

        //cout << "build tr... " <<endl ;

        unsigned int n_visit_mark = 0, curIdx = 0;
        visit_mark[timeStamp][n_visit_mark++] = uStart;
        visit[timeStamp][uStart] = true;
        hyperG[timeStamp][uStart].push_back(hyperId);
        //++spread[uStart];

        if (inDeg[timeStamp][uStart] == 0){
            // if uStart is a single node/has no in-neighbors, not in while loop
            curIdx++;
        }


        while (curIdx < n_visit_mark) // IC model
        {
            int i = visit_mark[timeStamp][curIdx++];
            // Traverse all in-degree neighbors and filter neighbors with timeStamp
            for (const auto& edge : T_graph[timeStamp]) {
                if (edge.ed == i){
                    int v = edge.st;  // in_neighbor id
                    double weight = edge.w;  // weight
                    double randDouble = sfmt_genrand_real1(&sfmtSeed);
                    if (visit[timeStamp][v])continue;
                    if (randDouble > weight)continue;

                    visit[timeStamp][v] = true;
                    visit_mark[timeStamp][n_visit_mark++] = v;
                    hyperG[timeStamp][v].push_back(hyperId);

                }
            }

        }

        for (unsigned int i = 0; i < n_visit_mark; i++) visit[timeStamp][visit_mark[timeStamp][i]] = false; //方便构造下个rr set

        hyperGT[timeStamp].push_back(vector<int>(visit_mark[timeStamp].begin(), visit_mark[timeStamp].begin() + n_visit_mark));


        return;
    }




    void buildhypergraphSCC(unsigned int R, const int timeStamp,
                            const unordered_map<int, unordered_map<int,vector<int>>>& FT_sccs,
                            const unordered_map<int, unordered_map<int, vector<pair<int, double>>>>& fin_all_scc_dags) // build TR sets of graph at timeStamp
    {

        if (R > INT_MAX){
            cout << "Error:R too large" << endl;
            exit(1);
        }
        //ASSERT(hyperGT.size() == 0);
        unsigned int size = hyperGT[timeStamp].size();
        int sccStart;  // SCC ID

        const auto& sccs = FT_sccs.at(timeStamp);
        const auto& scc_dag = fin_all_scc_dags.at(timeStamp);
        vector<int> sizes;
        sizes.reserve(sccs.size());
        for (const auto& scc_pair : sccs) {
            sizes.push_back(scc_pair.second.size());
        }

        // 前缀和 + 二分查找
        std::vector<double> prefix_sum(sizes.size());
        prefix_sum[0] = sizes[0];
        for (int i = 1; i < sizes.size(); ++i)
            prefix_sum[i] = prefix_sum[i - 1] + sizes[i];

        while (size < R){
            double r = sfmt_genrand_real1(&sfmtSeed) * prefix_sum.back();
            int sccStart = std::lower_bound(prefix_sum.begin(), prefix_sum.end(), r) - prefix_sum.begin();
            BuildHypergraphNodeSCC(sccStart, size++, timeStamp, sccs, scc_dag);

        }

    }

    void BuildHypergraphNodeSCC(const int sccStart, unsigned int hyperId, const int timeStamp,
                                const unordered_map<int, vector<int>>& sccs,
                                const unordered_map<int, vector<pair<int, double>>> & scc_dag)
    {

        //cout << "build RR-SCC on scc-dag, sampled scc id is: " << sccStart <<endl ;

        unsigned int n_visit_mark = 0, curIdx = 0;
        visit_mark[timeStamp][n_visit_mark++] = sccStart;
        visit[timeStamp][sccStart] = true;
        for (int u : sccs.at(sccStart)) {
            hyperG[timeStamp][u].push_back(hyperId);
        }
        //++spread[uStart];

        if (scc_dag.count(sccStart) == 0 || scc_dag.at(sccStart).empty()){
             //cout << "scc_dag size: 0" << endl;
            // if uStart is a single node/has no in-neighbors, not in while loop
            curIdx++;
        }


        while (curIdx < n_visit_mark) // IC model
        {
            int sccid = visit_mark[timeStamp][curIdx++];
            // Traverse all in-degree neighbors and filter neighbors with timeStamp
            auto it = scc_dag.find(sccid);
            if (it != scc_dag.end()) {
                for (const auto& edge : it->second) {
                    int scc_v = edge.first;  // in_neighbor scc id
                    double weight = edge.second;  // weight
                    double randDouble = sfmt_genrand_real1(&sfmtSeed);
                    if (visit[timeStamp][scc_v])continue;
                    //cout << "random: " <<randDouble << " weight: " << weight <<endl;
                    if (randDouble > weight)continue;

                    visit[timeStamp][scc_v] = true;
                    visit_mark[timeStamp][n_visit_mark++] = scc_v;
                    //hyperG[timeStamp][scc_v].push_back(hyperId);
                    //const auto& nodes = sccs[scc_v];
                    const auto& nodes = sccs.at(scc_v);
                    for (int u : nodes) {
                        hyperG[timeStamp][u].push_back(hyperId);
                    }
//                    for (int u : sccs[scc_v]) {
//                        hyperG[timeStamp][u].push_back(hyperId);
//                    }
                }
            }
        }

        for (unsigned int i = 0; i < n_visit_mark; i++) visit[timeStamp][visit_mark[timeStamp][i]] = false; //方便构造下个rr set

        unordered_set<int> node_set;
        for (int i = 0; i < n_visit_mark; ++i) {
            int scc_id = visit_mark[timeStamp][i];
            const auto& nodes = sccs.at(scc_id);
            for (int u : nodes) {
                node_set.insert(u);
            }
        }
        //cout << "node_set.size: " << node_set.size() << endl;
        hyperGT[timeStamp].push_back(vector<int>(node_set.begin(), node_set.end()));

        return;
    }





//// ================================== Basic Function =======================================

// 得到所有时刻T_valid
    void get_T_valid(vector<int> &T_valid){
        for (int i = 0; i < T_graph.size(); i ++) {
            if (!T_graph[i].empty()) {

                T_valid.push_back(i);
            }
        }
    }

    void getT_graph_TargetSet() // 得到每个时刻图下的目标节点集合
    {
        cout << "Get Target node set of each snapshot graph ..." << endl;
        // 遍历每个时刻
        for (int t_snap = 0; t_snap < T_graph.size(); ++t_snap) {
            set<int> targetsAtTime;  // 存储当前时刻的目标节点
            // 遍历当前时刻下的每个节点
            if (!T_graph[t_snap].empty()){
                for (const auto& node_edge : T_graph[t_snap]) {
                    int u = node_edge.st;  // 获取节点
                    int v = node_edge.ed;

                    targetsAtTime.insert(u);
                    targetsAtTime.insert(v);
                }
                // 将当前时刻的目标节点集合加入到targetsOfT_graph中
                targetsOfT_graph.insert(make_pair(t_snap, targetsAtTime));
            }
            //cout << "time " << t_snap << " numner is: " << targetsAtTime.size()<<endl;
        }
    }



//// ========= Persistence sequence Computation =====================
// 辅助函数，找到所有连续的子序列，且尾 - 首 +1>= timespan
    void findValidSequences(const vector<int> seq, int theta, vector<pair<int, int>>& results) {
        int n = seq.size();
        if (n == 0) return;

        int start = 0;  // 当前连续子序列的起始位置
        for (int i = 1; i <= n; i++) {
            // 如果到达了序列末尾或者当前元素不再连续
            if (i == n || seq[i] != seq[i - 1] + 1) {
                // 计算当前连续子序列的长度
                int length = i - start;
                //cout << "Checking subsequence: Start: " << seq[start] << ", End: " << seq[i - 1] << ", Length: " << length << endl;
                //cout << "length: "<<length <<" timespan: "<< theta<<endl;
                if (length >= theta) {
                    // 满足条件，存储首尾
                    //cout << "Adding to results: Start: " << seq[start] << ", End: " << seq[i - 1] << ", Length: " << length << endl;
                    results.push_back(make_pair(seq[start], seq[i - 1]));
                }
                // 更新起始位置
                start = i;
            }
        }
        /*
        for (const auto& r : results) {
            cout << "Start: " << r.first << ", End: " << r.second << endl;
        }*/
    }

// 找到所有窗口内连续的子序列
    void findValidSequencesWin(const vector<int> seq, const int theta, const int W, vector<pair<int, int>>& results) {
        int n = seq.size();
        if (n == 0) return;

        int interval_ts=-1, interval_te=-1;
        vector<int> prefixSum(Tmax + 1, 0);  // 前缀和

        //// 填充前缀和数组
        int lastSum = 0;
        int lastTime = Tmin - 1;
        for (int t : seq) {
            fill(prefixSum.begin() + lastTime + 1, prefixSum.begin() + t + 1, lastSum);
            prefixSum[t] = lastSum + 1; // 只在 isAlpha[t] == true 处增加
            lastSum = prefixSum[t];
            lastTime = t;
        }
        fill(prefixSum.begin() + lastTime + 1, prefixSum.end(), lastSum);
//                    std::cout << "Prefix Sum Array: ";
//                    for (int i = Tmin; i <= Tmax; ++i) {
//                        std::cout << prefixSum[i] << " ";
//                    }
//                    std::cout << std::endl;

        // 使用滑动窗口法，遍历每个可能的窗口
        int Tmin_time=max(0, seq[theta - 1]);;
        int Tmax_time=seq[seq.size()-1];
        //cout << Tmin_node << ", " << Tmax_node<< endl;
        for (int i = Tmin_time - W+1; i <= Tmax_time + W+1; i++) {
            if (i<0) continue;
            if (i + W-1>Tmax_time) continue;
            // 计算窗口 [i, i+W-1] 内的有效时刻数
            int count = prefixSum[i + W-1] - prefixSum[i-1];
            //cout << "start: "<< i << " end: " << i + W-1 << " alpha num: " << count <<endl;//getchar();
            if (count >= theta) {
                //cout << "start: "<< i << " end: " << i + W-1 << " alpha num: " << count <<endl; //getchar();
                if (interval_ts == -1 && interval_te == -1){ //第一个有效的时间窗口，还未开始连续窗口
                    interval_ts = i;
                    interval_te = i + W-1;
                }
                else if (interval_te+1 == i + W -1){ //找到连续窗口，ts不变，te+1
                    interval_te = i + W-1;
                }
            }
            else { // 遇到不满足theta的窗口，记录目前的最大时间段
                if (interval_ts != -1 && interval_te != -1) {
                    results.push_back(make_pair(interval_ts, interval_te));
                }

                interval_ts = -1;
                interval_te = -1;
            }
        }
        // 处理最后一个有效时间段（如果存在）
        if (interval_ts != -1 && interval_te != -1) {
            results.push_back(make_pair(interval_ts, interval_te));
        }
        /*
        for (const auto& r : results) {
            cout << "Start: " << r.first << ", End: " << r.second << endl;
        }*/
    }

// 判断pers_t[node]中的时刻是否存在满足条件的时刻窗口
    bool hasThetaInRange(const vector<int>& T, int W, int theta) {
        int left = 0;
        int n = T.size();

        for (int right = 0; right < n; ++right) {
            while (T[right] - T[left] > W) {
                left++;
            }
            if (right - left + 1 >= theta) {
                return true;  // 找到满足条件的窗口
            }
        }
        return false;
    }



    //// ========= Combination method =====================
// 递归辅助函数，用来生成组合
    void generateCombinations(const std::vector<int>& non_target, int k, int start, std::vector<int>& current_combination, std::vector<std::vector<int>>& seedsets_candidate) {
        if (current_combination.size() == k) {
            // 将当前组合加入候选种子集
            seedsets_candidate.push_back(vector<int>(current_combination.begin(), current_combination.end()));
            return;
        }

        // 从 start 位置开始选择元素
        for (int i = start; i < non_target.size(); ++i) {
            current_combination.push_back(non_target[i]);
            generateCombinations(non_target, k, i + 1, current_combination, seedsets_candidate);
            current_combination.pop_back(); // 回溯
        }
    }



//// ========= ReverseGreedy functions =====================

    vector<int> select_topk_nodes_by_degree(const vector<edge>& T_graph_i, int k) {
        unordered_map<int, int> out_deg;
        for (const auto& e : T_graph_i) {
            out_deg[e.st]++;
        }

        vector<pair<int, int>> nodes(out_deg.begin(), out_deg.end());
        sort(nodes.begin(), nodes.end(), [](auto& a, auto& b) {
            return a.second > b.second;
        });

        vector<int> topk;
        for (int i = 0; i < min(k, (int)nodes.size()); ++i) {
            topk.push_back(nodes[i].first);
        }
        return topk;
    }

    vector<int> greedy_select_topk_nodes_by_rr(int t_snap, int k, int R) {
        vector<bool> rr_covered(R, false);  // 标记哪些 RR 已被覆盖
        unordered_map<int, unordered_set<int>> node_rr_map; // node -> RR ids
        unordered_map<int, int> node_gain; // 当前节点能覆盖的 RR 数

        // 初始化 node -> RR ids 映射（从 hyperG 中来）
        for (int node : nodeset) {
            auto node_it = hyperG[t_snap].find(node);
            if (node_it == hyperG[t_snap].end()) continue;
            for (int rr_id : node_it->second) {
                node_rr_map[node].insert(rr_id);
                node_gain[node]++;
            }
        }


        vector<int> topk;
        for (int i = 0; i < k; ++i) {
            // 1. 选择当前边际增益最大的节点
            int best_node = -1, best_gain = -1;
            for (auto it = node_gain.begin(); it != node_gain.end(); ++it) {
                int node = it->first;
                double gain = it->second;
                if (gain > best_gain) {
                    best_node = node;
                    best_gain = gain;
                }
            }
            if (best_node == -1 || best_gain == 0) break;  // 提前结束

            //cout << "current best_node: " << best_node << " gain: " << best_gain<<endl;
            topk.push_back(best_node);

            // 2. 更新 RR 覆盖情况，并更新其他节点的 gain
            for (int rr_id : node_rr_map[best_node]) {
                if (!rr_covered[rr_id]) {
                    rr_covered[rr_id] = true;
                    // 对所有包含该 RR 的节点，其 gain 减一
                    for (int node : hyperGT[t_snap][rr_id]) {
                        if (node != best_node && node_gain.count(node)) {
                            node_gain[node]--;
                        }
                    }
                }
            }

            // 从 node_gain 中移除已选节点
            node_gain.erase(best_node);
        }
//        for (int i : topk) {
//                cout << i << " ";  // 输出每个值
//        }cout <<"]"<<endl <<endl;

        return topk;
    }


//// ========= LazyReplace functions =====================



//// ========= Other functions =====================


    int find_best_initset(vector<int> & fin_can_seedset, vector<int> &fin_remain_nodes,
                          const int k, const int R, const double alpha, const int theta){

        // inf top-k
        vector<int> can_seedset1;
        vector<int> remain_nodes1;
        unordered_map<int, double> singnode_inf;
        for (int node : nodeset) {
            for (int t_snap : T_valid) {
                auto snap_it = hyperG.find(t_snap);
                if (snap_it == hyperG.end()) continue;

                auto node_it = snap_it->second.find(node);
                if (node_it == snap_it->second.end()) continue;

                int coverTR = node_it->second.size();
                double inf_t = ((double)coverTR / R) * (double)size_n;
                singnode_inf[node] += inf_t;
            }
        }
        vector<int> extra_nodes;
        // 直接按照avg_inf大小对节点排序, 选前top-k个节点作为初始can_seedset
        for (int i = 0; i < size_n; ++i) {
            extra_nodes.push_back(i);

        }
        sort(extra_nodes.begin(), extra_nodes.end(),
             [&singnode_inf](int a, int b) {  // 这里改为 int 比较
                 return singnode_inf[a] > singnode_inf[b];
             });
        for (int i = 0; i < extra_nodes.size(); ++i) {
            if (i < k) can_seedset1.push_back(extra_nodes[i]);
            //else if (i < int(alpha*2 * size_n)+k) remain_nodes1.push_back(extra_nodes[i]);
            else remain_nodes1.push_back(extra_nodes[i]);
        }
        set<int> curr_perst;
        vector<pair<int, int>> persistent_inf_interval;
        int inf_persistence1 = 0;
        int coverTR;
        double inf_t;

        for (int i = 0; i < T_valid.size(); i++){
            int t_snap = T_valid[i];
            //double inf_t = tg.GetInfTgraph_Set_MC(R, can_seedset, t_snap);
            vector<bool> rr_bitmap(R, false);
            for (int node : can_seedset1) {
                for (int idx : hyperG[t_snap][node]) {
                    rr_bitmap[idx] = true;
                }
            }
            coverTR = std::count(rr_bitmap.begin(), rr_bitmap.end(), true);

            inf_t =(( double) coverTR / R) * (double) size_n;
            //cout << "inf_t: " << inf_t << " coverTR: " << coverTR <<endl;

            if (inf_t >= alpha * size_n){
                curr_perst.insert(t_snap);
            }
        }

        cout <<"pers_t.size is " << curr_perst.size() <<endl;

        if (curr_perst.size() >= theta) {
            vector<int> seq(curr_perst.begin(), curr_perst.end());
            findValidSequences(seq, theta, persistent_inf_interval);
            cout << "number of valid interval："<< persistent_inf_interval.size() <<endl;
            for (const auto& p : persistent_inf_interval) {
                int length = p.second - p.first + 1;  // 计算长度
                inf_persistence1 += length;     // 累加该节点所有序列的长度
                //cout << "Sequence (head: " << p.first << ", tail: " << p.second << ") "<<endl;
            }
            cout << "inf with the total sequence length: " << inf_persistence1 << endl;
        }

        // freq top-k
        vector<int> can_seedset2;
        vector<int> remain_nodes2;
        unordered_map<int, vector<int>> snapshot_seedsets; // 每个 snapshot 对应的 Si
        cout << "selected each snapshot graph top-k user by greedy" <<endl;
        for (int i: T_valid) {
            snapshot_seedsets[i] = greedy_select_topk_nodes_by_rr(i, k, R);
        }

        vector<int> node_freq (size_n, 0);
        for (int t : T_valid) {
            for (int node : snapshot_seedsets[t]) {
                node_freq[node]++;
            }
        }

        extra_nodes.clear();
        // 直接按照频率大小对节点排序, 选前top-k个节点作为初始can_seedset
        for (int i = 0; i < size_n; ++i) {
            extra_nodes.push_back(i);

        }
        sort(extra_nodes.begin(), extra_nodes.end(),
             [&node_freq](int a, int b) {  // 这里改为 int 比较
                 return node_freq[a] > node_freq[b];
             });
        for (int i = 0; i < extra_nodes.size(); ++i) {
            if (i < k) can_seedset2.push_back(extra_nodes[i]);
            //else if (i < int(alpha *2 * size_n)+k) remain_nodes2.push_back(extra_nodes[i]);
            else remain_nodes2.push_back(extra_nodes[i]);
        }

        curr_perst.clear();
        persistent_inf_interval.clear();
        int inf_persistence2 = 0;

        //// compute fist-joint influence of can_seedset
        for (int i = 0; i < T_valid.size(); i++){
            int t_snap = T_valid[i];
            //double inf_t = tg.GetInfTgraph_Set_MC(R, can_seedset, t_snap);
            vector<bool> rr_bitmap(R, false);
            for (int node : can_seedset2) {
                for (int idx : hyperG[t_snap][node]) {
                    rr_bitmap[idx] = true;
                }
            }
            coverTR = std::count(rr_bitmap.begin(), rr_bitmap.end(), true);

            inf_t =(( double) coverTR / R) * (double) size_n;
            //cout << "inf_t: " << inf_t << " coverTR: " << coverTR <<endl;

            if (inf_t >= alpha * size_n){
                curr_perst.insert(t_snap);
            }
        }

        cout <<"pers_t.size is " << curr_perst.size() <<endl;
//        cout << "time stamps satisfying constraints:[ " ;
//        for (int val : curr_perst) {
//            cout << val << " ";  // 输出每个w值
//        }cout <<"]"<<endl;
        if (curr_perst.size() >= theta) {
            vector<int> seq(curr_perst.begin(), curr_perst.end());
            findValidSequences(seq, theta, persistent_inf_interval);
            cout << "number of valid interval："<< persistent_inf_interval.size() <<endl;
            for (const auto& p : persistent_inf_interval) {
                int length = p.second - p.first + 1;  // 计算长度
                inf_persistence2 += length;     // 累加该节点所有序列的长度
//                cout << "Sequence (head: " << p.first << ", tail: " << p.second << ") "<<endl;
            }
            cout << "freq with the total sequence length: " << inf_persistence2 << endl;
        }

        if (inf_persistence2 >= inf_persistence1){
            fin_can_seedset = can_seedset2;
            fin_remain_nodes = remain_nodes2;
            return inf_persistence2;
        }
        else{
            fin_can_seedset = can_seedset1;
            fin_remain_nodes = remain_nodes1;
            return inf_persistence1;
        }

    }


    int find_best_initset_win(vector<int> & fin_can_seedset, vector<int> &fin_remain_nodes,
                          const int k, const int R, const double alpha, const int theta, const int W){

        // inf top-k
        vector<int> can_seedset1;
        vector<int> remain_nodes1;
        unordered_map<int, double> singnode_inf;
        for (int node : nodeset) {
            for (int t_snap : T_valid) {
                auto snap_it = hyperG.find(t_snap);
                if (snap_it == hyperG.end()) continue;

                auto node_it = snap_it->second.find(node);
                if (node_it == snap_it->second.end()) continue;

                int coverTR = node_it->second.size();
                double inf_t = ((double)coverTR / R) * (double)size_n;
                singnode_inf[node] += inf_t;
            }
        }
        vector<int> extra_nodes;
        // 直接按照avg_inf大小对节点排序, 选前top-k个节点作为初始can_seedset
        for (int i = 0; i < size_n; ++i) {
            extra_nodes.push_back(i);

        }
        sort(extra_nodes.begin(), extra_nodes.end(),
             [&singnode_inf](int a, int b) {  // 这里改为 int 比较
                 return singnode_inf[a] > singnode_inf[b];
             });
        for (int i = 0; i < extra_nodes.size(); ++i) {
            if (i < k) can_seedset1.push_back(extra_nodes[i]);
            //else if (i < int(alpha * size_n)+k) remain_nodes1.push_back(extra_nodes[i]);
            else remain_nodes1.push_back(extra_nodes[i]);
        }
        set<int> curr_perst;
        vector<pair<int, int>> persistent_inf_interval;
        int inf_persistence1 = 0;
        int coverTR;
        double inf_t;

        for (int i = 0; i < T_valid.size(); i++){
            int t_snap = T_valid[i];
            //double inf_t = tg.GetInfTgraph_Set_MC(R, can_seedset, t_snap);
            vector<bool> rr_bitmap(R, false);
            for (int node : can_seedset1) {
                for (int idx : hyperG[t_snap][node]) {
                    rr_bitmap[idx] = true;
                }
            }
            coverTR = std::count(rr_bitmap.begin(), rr_bitmap.end(), true);

            inf_t =(( double) coverTR / R) * (double) size_n;
            //cout << "inf_t: " << inf_t << " coverTR: " << coverTR <<endl;

            if (inf_t >= alpha * size_n){
                curr_perst.insert(t_snap);
            }
        }

        cout <<"pers_t.size is " << curr_perst.size() <<endl;

        if (curr_perst.size() >= theta) {
            vector<int> seq(curr_perst.begin(), curr_perst.end());
            findValidSequencesWin(seq, theta, W, persistent_inf_interval);
            cout << "number of valid interval："<< persistent_inf_interval.size() <<endl;
            for (const auto& p : persistent_inf_interval) {
                int length = p.second - p.first + 1;  // 计算长度
                inf_persistence1 += length;     // 累加该节点所有序列的长度
                //cout << "Sequence (head: " << p.first << ", tail: " << p.second << ") "<<endl;
            }
            cout << "inf with the total sequence length: " << inf_persistence1 << endl;
        }

        // freq top-k
        vector<int> can_seedset2;
        vector<int> remain_nodes2;
        unordered_map<int, vector<int>> snapshot_seedsets; // 每个 snapshot 对应的 Si
        cout << "selected each snapshot graph top-k user by greedy" <<endl;
        for (int i: T_valid) {
            snapshot_seedsets[i] = greedy_select_topk_nodes_by_rr(i, k, R);
        }

        vector<int> node_freq (size_n, 0);
        for (int t : T_valid) {
            for (int node : snapshot_seedsets[t]) {
                node_freq[node]++;
            }
        }

        extra_nodes.clear();
        // 直接按照频率大小对节点排序, 选前top-k个节点作为初始can_seedset
        for (int i = 0; i < size_n; ++i) {
            extra_nodes.push_back(i);

        }
        sort(extra_nodes.begin(), extra_nodes.end(),
             [&node_freq](int a, int b) {  // 这里改为 int 比较
                 return node_freq[a] > node_freq[b];
             });
        for (int i = 0; i < extra_nodes.size(); ++i) {
            if (i < k) can_seedset2.push_back(extra_nodes[i]);
            //else if (i < int(alpha * size_n)+k) remain_nodes2.push_back(extra_nodes[i]);
            else remain_nodes2.push_back(extra_nodes[i]);
        }

        curr_perst.clear();
        persistent_inf_interval.clear();
        int inf_persistence2 = 0;

        //// compute fist-joint influence of can_seedset
        for (int i = 0; i < T_valid.size(); i++){
            int t_snap = T_valid[i];
            //double inf_t = tg.GetInfTgraph_Set_MC(R, can_seedset, t_snap);
            vector<bool> rr_bitmap(R, false);
            for (int node : can_seedset2) {
                for (int idx : hyperG[t_snap][node]) {
                    rr_bitmap[idx] = true;
                }
            }
            coverTR = std::count(rr_bitmap.begin(), rr_bitmap.end(), true);

            inf_t =(( double) coverTR / R) * (double) size_n;
            //cout << "inf_t: " << inf_t << " coverTR: " << coverTR <<endl;

            if (inf_t >= alpha * size_n){
                curr_perst.insert(t_snap);
            }
        }

        cout <<"pers_t.size is " << curr_perst.size() <<endl;
//        cout << "time stamps satisfying constraints:[ " ;
//        for (int val : curr_perst) {
//            cout << val << " ";  // 输出每个w值
//        }cout <<"]"<<endl;
        if (curr_perst.size() >= theta) {
            vector<int> seq(curr_perst.begin(), curr_perst.end());
            findValidSequencesWin(seq, theta, W, persistent_inf_interval);
            cout << "number of valid interval："<< persistent_inf_interval.size() <<endl;
            for (const auto& p : persistent_inf_interval) {
                int length = p.second - p.first + 1;  // 计算长度
                inf_persistence2 += length;     // 累加该节点所有序列的长度
//                cout << "Sequence (head: " << p.first << ", tail: " << p.second << ") "<<endl;
            }
            cout << "freq with the total sequence length: " << inf_persistence2 << endl;
        }

        if (inf_persistence2 >= inf_persistence1){
            fin_can_seedset = can_seedset2;
            fin_remain_nodes = remain_nodes2;
            return inf_persistence2;
        }
        else{
            fin_can_seedset = can_seedset1;
            fin_remain_nodes = remain_nodes1;
            return inf_persistence1;
        }

    }


//// ========= Output information =====================

    void outTimeStamps_NodeSatisfyAlpha (vector<int> pers_t_node){
        cout << "time stamps satisfying constraints:[ " ;
        for (int val : pers_t_node) {
            cout << val << " ";  // 输出每个值
        }cout <<"]"<<endl;
    }

    void seedset_validation (vector<int> final_seedset, vector<int> T, int R, double alpha, double theta){
        set<int> val_alpha_perst;

        std::ofstream outfile("../../casestudy/output-inf-EC-t3.txt");
        for (int i :T) {
            vector<bool> rr_bitmap(R, false);
            for (int node : final_seedset) {
                for (int idx : hyperG[i][node]) {
                    rr_bitmap[idx] = true;
                }
            }
            int coverTR = std::count(rr_bitmap.begin(), rr_bitmap.end(), true);

            double inf_t =(( double) coverTR / R) * (double) size_n;
            //cout << "inf_t: " << inf_t << " coverTR: " << coverTR <<endl;

            if (inf_t >= alpha * size_n)
                val_alpha_perst.insert(i);

            double ratio = inf_t / size_n;
            outfile << i-Tmin << " " << std::fixed << std::setprecision(4) << ratio << std::endl;

        }
        outfile.close();

        cout << "[val alpha] total number of time stamps: " << val_alpha_perst.size()<<endl;

        set<int> val_theta_perst;
        int val_init_perstence = 0;
        vector<pair<int, int>> val_persistent_inf_interval;
        if (val_alpha_perst.size() >= theta) {
            vector<int> seq(val_alpha_perst.begin(), val_alpha_perst.end());
            // 查找所有连续子序列且满足尾 - 首 +1 >= timespan
            findValidSequences(seq, theta, val_persistent_inf_interval);
            // 计算每个有效序列的长度并累加
            cout << "number of validation valid interval："<< val_persistent_inf_interval.size() <<endl;

            std::ofstream outseq("../../casestudy/seq-inf-EC-t3.txt");
            for (const auto& p : val_persistent_inf_interval) {
                int length = p.second - p.first +1;
                val_init_perstence += length;
                for (int t = p.first; t <=p.second; t++){
                    val_theta_perst.insert(t);
                }
                outseq << p.first-Tmin << " " << p.second-Tmin << std::endl;
                cout << "Sequence (head: " << p.first << ", tail: " << p.second << ") "<<endl;
            }
            outseq.close();

        }
        cout << "[val theta] total number of time stamps: " << val_theta_perst.size()<<endl;
        cout << "validate total inf_perstence: " << val_init_perstence<<endl;

    }

    void seedset_validationWin (vector<int> final_seedset, vector<int> T, int R, double alpha, double theta, int W){
        set<int> val_alpha_perst;

        std::ofstream outfile("../../casestudy/outputW-EC-t5-W7.txt");
        for (int i :T) {
            vector<bool> rr_bitmap(R, false);
            for (int node : final_seedset) {
                for (int idx : hyperG[i][node]) {
                    rr_bitmap[idx] = true;
                }
            }
            int coverTR = std::count(rr_bitmap.begin(), rr_bitmap.end(), true);

            double inf_t =(( double) coverTR / R) * (double) size_n;
            //cout << "inf_t: " << inf_t << " coverTR: " << coverTR <<endl;

            if (inf_t >= alpha * size_n)
                val_alpha_perst.insert(i);

            double ratio = inf_t / size_n;
            outfile << i-Tmin << " " << std::fixed << std::setprecision(4) << ratio << std::endl;

        }
        outfile.close();

        cout << "[val alpha] total number of time stamps: " << val_alpha_perst.size()<<endl;

        set<int> val_theta_perst;
        int val_init_perstence = 0;
        vector<pair<int, int>> val_persistent_inf_interval;
        if (val_alpha_perst.size() >= theta) {
            vector<int> seq(val_alpha_perst.begin(), val_alpha_perst.end());
            findValidSequencesWin(seq, theta, W, val_persistent_inf_interval);
            cout << "number of validation valid interval："<< val_persistent_inf_interval.size() <<endl;

            std::ofstream outseq("../../casestudy/seqW-EC-t5-W7.txt");
            for (const auto& p : val_persistent_inf_interval) {
                int length = p.second - p.first +1;
                val_init_perstence += length;
                for (int t = p.first; t <=p.second; t++){
                    val_theta_perst.insert(t);
                }
                outseq << p.first-Tmin << " " << p.second-Tmin << std::endl;
                cout << "Sequence (head: " << p.first << ", tail: " << p.second << ") "<<endl;
            }
            outseq.close();

        }
        cout << "[val theta] total number of time stamps: " << val_theta_perst.size()<<endl;
        cout << "validate total inf_perstence: " << val_init_perstence<<endl;
    }








};



















