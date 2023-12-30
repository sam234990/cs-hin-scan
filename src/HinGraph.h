#ifndef _HINGRAPH_H_
#define _HINGRAPH_H_

#include "utils.h"
#include "timer.h"
#include "multiLevelQueue.h"
#include "MuTree.h"
#include <vector>
#include <string>
#include <bitset>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <queue>
#include <cmath>
#include <thread>
#include <functional>  // 添加这个头文件

struct Vertex_neighbor
{
    unordered_map<int, vector<int>> d_neighbor_;
    unordered_map<int, int> num_d_neighbor;
    vector<int> type_order_;
};

struct Nei_similarity
{
    Nei_similarity(int nei_id, float id_sim) : neighbor_i(nei_id)
    {
        sim_vec.push_back(id_sim);
    }
    Nei_similarity() : neighbor_i(0), domin_rank(0) {}

    int neighbor_i;
    vector<float> sim_vec;
    int domin_rank;
};

struct Query_nei_similarity
{
    Query_nei_similarity(int nei_id, float id_sim) : neighbor_i(nei_id), similarity(id_sim) {}

    int neighbor_i;
    float similarity;
};

const int SP_SIZE = 100;

struct k_threshold
{
    vector<vector<float>> thres_vecs;
    vector<vector<float>> corner_points; // contain the thres of dim0 and dim1, sort by the descending order of dim1
    vector<bool> used_neighbor;
    bitset<SP_SIZE + 1> fix_thres_dim1;
    bitset<SP_SIZE + 1> unsat_dim1;
    bitset<SP_SIZE + 2> one_dim;
    int max_dim0;
};

struct k_homo_adj_node
{
    k_homo_adj_node() : neighbor_i(0), re_type_sim(0.0) {}
    k_homo_adj_node(int id, float sim) : neighbor_i(id), re_type_sim(sim) {}
    int neighbor_i;
    float re_type_sim;
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
    vector<vector<int>> qn_adj_List; // qn
    vector<int> pa;               // pa and rank are used for the disjoint-set data structure
    vector<int> p_rank_;

    // index variables
    int k_max;
    vector<int> index_type_order;
    vector<float> index_order_epsilon;
    vector<bool> last_vertex;
    vector<bool> community_vertex;
    vector<Vertex_neighbor> dn_adj_List; // dn
    vector<vector<Nei_similarity>> h_sim;
    vector<k_threshold> node_k_thres;
    unordered_map<int, vector<int>> t_hop_visited;
    unordered_map<int, vector<int>> empty_dn_set;
    vector<vector<k_homo_adj_node>> k_homo_graph;

    // imporve_index variables
    int all_community_num;
    int has_community;

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
    vector<bool> res_community;
    multiLevelQueue mlq;

    void load_query_file(string query_file);
    void initialize_query_();
    void reinitialize_query_();
    void baseline_query_();
    void online_query_scan();
    void scan_check_cluster_core(int u);

    // void improved_query_();

    // void estimate_best_order(vector<int> &order_type);

    void print_result(bool print_all, long use_time);

    // cs-search
    void search_k_strata(int u);
    void search_cand_sn(int u);
    void check_sn(int u, queue<int> &cs_queue, queue<int> &delete_q);
    void core_decomposition(queue<int> &delete_q);

    // index construct
    void initial_construct_index();
    void search_d_neighbor();
    void save_d_n();
    void load_d_n();

    void multi_func(int i);
    void compute_all_similarity();
    int compute_one_type_qn_similarity(int i, int type_i, vector<Query_nei_similarity> &type_i_qn_similarity);
    void concat_one_type_qn(const vector<int> &dn_i_type_i, vector<int> &mVector, bool query_type);
    void intersection_neisim(vector<Nei_similarity> &nei_sim, const vector<Query_nei_similarity> &type_i_qn_similarity);
    void compute_domin_rank(vector<Nei_similarity> &qn_sim);
    void save_all_similarity();
    void load_all_similarity();

    void compute_k_threshold(int start_k);
    void compute_connect_k_core(int k, const vector<float> &fix_type, int re_type);
    bool search_and_add_threshold(int k, const vector<float> &fix_type, float re_type_threshold);
    void save_k_thres_vec(int k);
    void load_k_thres_vec(int k);

    // imporve index construct
    void improve_k_thres(int start_k);
    void skyline_highD(int k, vector<float> &cons, int cur_d);
    void skyline3D(int k, vector<float> &cons);
    void skyline2D(int k, const vector<float> cons);
    // void update_concer_point(const vector<float> &cons, k_threshold &thres_corner, int type_i);

    float constraint_one_dim(int k, const vector<float> cons, int type_i, int fix_vertex_thres_);
    bool compute_one_dim_max(int k, float re_type_threshold, const vector<float> cons,
                             int type_i, const vector<int> visit, vector<bool> fix_vertex);
    bool bfs_community(int start_i, vector<int> &visit, const vector<bool> fix_vertex,
                       int community_num, bitset<SP_SIZE + 1> &c_dim1);

    // MuTree index
    void build_tree_index(int start_k);
    void build_one_type_tree(int k, int type_j, int cons);

    // index query
    void index_query_();
    void query_index_scan();
    bool index_judge_core(int i, int k);
    void select_query_node();

    // similarity compute
    bool check_struc_sim(int a, int b);
    int get_vertex_type(int vertex_id);
    double calculateJaccardSimilarity(const vector<int> &vec1, const vector<int> &vec2);
    bool judgeJacSim(const vector<int> &vec1, const vector<int> &vec2, double type_i_epsilon);

    void my_union(int u, int v);
    int find_root(int u);

    // debug check
    void check_empty_set();
    void count_core_vertices();

    void dfs(int cur_i, int step, int types, vector<long> &res);

public:
    // index variables
    vector<MuTree> index_tree;

    HinGraph(string data_dir);
    ~HinGraph();

    void load_graph();
    void output_result(string output);
    void cs_hin_scan(string query_file, string mode);
    void construct_index(string query_file, string option, int start_k, int scale);
    void find_meta(int type);

    int max_scale_num, max_scale_id;
    int scale_;
    string option_;
};

#endif