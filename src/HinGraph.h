#ifndef _HINGRAPH_H_
#define _HINGRAPH_H_

#include "utils.h"
#include "timer.h"
#include "MuTree.h"
#include "ONIndex.h"
#include <vector>
#include <string>
#include <vector>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <queue>
#include <cmath>

struct Vertex_neighbor
{
    unordered_map<int, vector<int>> d_neighbor_;
    unordered_map<int, int> num_d_neighbor;
};

struct Query_neighbor_order
{
    // vector<Query_nei_similarity> type_nei_order_;
    // vector<int> start_offset_;
    vector<Query_nei_similarity> type_nei_order_;
    // int type_nei_num;
    // vector<vector<Query_nei_similarity>> type_nei_order_;
    // vector<int> type_nei_num;
    // unordered_map<int, vector<Query_nei_similarity>> type_nei_order_;
};

class HinGraph
{
private:
    string data_file_name;
    string data_dir_;
    string data_index_dir_;
    string query_file_name;
    int mode_query;

    // hin graph
    int m, n, n_types;
    vector<vector<int>> hin_schema_adjacencyMatrix;
    vector<vector<int>> hin_schema_edge_cnt;
    vector<int> vertex_start_map_;
    int *vertex_offset_;
    Vertex_type *edges_;

    // query parameters
    int p_mu, p_d, p_query_type, query_i;
    int query_node_num;
    vector<int> query_node_list;
    map<int, double> type_epsilon;
    bool unit_epsilon, random_query;
    double unit_epsilon_value;

    // scan variables
    int query_type_offset_; // vertex_id = i + query_type_offset_, i is the subscript in adjacencyList
    int num_query_type_;
    map<int, int> distance_;
    vector<bool> cand_core_;
    vector<vector<int>> qn_adj_List;     // qn
    vector<Vertex_neighbor> dn_adj_List; // dn
    unordered_map<int, vector<int>> t_hop_visited;

    // pscan variables
    vector<bool> visited_qtv_;
    vector<int> similar_degree;
    vector<int> effective_degree;

    // cs variable
    vector<vector<int>> cand_sn_list; // candidate structure neighbor
    vector<vector<int>> non_sn_list;  // qn
    int cand_gen_type;
    vector<bool> ks_visit;
    vector<bool> is_in_community;

    void load_query_file(string query_file);
    void initialize_query_();
    void reinitialize_query_();
    void baseline_query_();
    void improved_query_();

    void estimate_best_order(vector<int> &order_type);

    void print_result(bool print_all);

    // cs-search
    void search_k_strata(int u);
    void search_cand_sn(int u);
    void cs_check_cluster_core(int u, queue<int> &cs_node_q);
    void explore_community(queue<int> &cs_queue);

    void cs_check_cluster_core_order(int u, queue<int> &cs_node_q, vector<int> order_type);

    // similarity compute
    bool check_struc_sim(int a, int b);
    int get_vertex_type(int vertex_id);
    double calculateJaccardSimilarity(const vector<int> &vec1, const vector<int> &vec2);
    bool judgeJacSim(const vector<int> &vec1, const vector<int> &vec2, double type_i_epsilon);
    bool check_one_type_sim(int a, int b, int type);
    bool check_struc_sim_with_order(int a, int b, const vector<int> order_type);

    // debug check
    void check_neighbor(int i);
    void check_two_hop_visited();
    void check_empty_set();
    void count_core_vertices();

public:
    // index variables
    vector<Query_neighbor_order> nei_orders_; // neighbor order with type
    ONIndex on_index_;
    map<double, int> unique_sim;
    vector<MuTree> index_tree;

    HinGraph(string data_dir);
    ~HinGraph();

    void load_graph();
    void output_result(string output);
    void cs_hin_scan(string query_file, string mode);
};

#endif