#pragma once

// #include "HinGraph.h"
#include <vector>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <algorithm> // For std::remove_if
#include <omp.h>     // 包含 OpenMP 头文件

using namespace std;

struct MetaPath
{
    vector<int> vertex;
    vector<int> edge;
    int pathLen;
};

struct pathcnt
{
    unordered_map<int, int> ins_path_cnt; // nei_offset -> pathcnt
};

struct path_hetesim
{
    unordered_map<int, double> hetesim; // neiid->hetesim
};

struct nei_list
{
    unordered_map<int, vector<int>> vertex_type_nei; // vertex_type -> [nei_id]
    unordered_map<int, vector<int>> edge_type_nei;   // edge_type -> [nei_id]
};

class HinGraph;

class PathSim
{
public:
    int path_num;
    std::vector<MetaPath> MetaPathVec;
    int query_num_;
    vector<vector<pathcnt>> query_path_cnt;
    vector<bool> path_search_finish_;
    vector<vector<int>> homo_graph;
    vector<int> homo_degree;
    vector<bool> p_induce;
    vector<vector<int>> p_g;

    // HeteSim
    vector<vector<path_hetesim>> hete_path_cnt;                    // [vertex_id][mp_id]
    vector<nei_list> vertex_nei_type_split; // [vertex_id][tycnt:[nei id]]
    vector<bool> OnlyVertexType;

    PathSim(/* args */);
    ~PathSim();
    void initial_metapaths(vector<string> mps);
    void initial_query_vertex(int q_num);
    void search(const HinGraph &graph, int query_i);
    double compute_avg_pathsim(const HinGraph &graph, int i, int j);
    bool judge_pathsim(int i, int j, const vector<float> sim_threshold);

    void generate_cand_nei(const HinGraph &graph, int query_i, vector<int> &cand_nei);

    void trans_homo_graph(const HinGraph &graph, string meta_path, string save_path);
    vector<int> p_induced_graph(const HinGraph &graph, int i);

    void initial_hetesim(const HinGraph &graph, int vertex_num, int query_num);
    bool judge_hetesim(int id1, int id2, const vector<float> sim_threshold);
    double compute_hetesim(int id1, int id2, int mp_id, int step);
    // void judge_hetesim()
};