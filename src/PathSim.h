#pragma once

// #include "HinGraph.h"
#include <vector>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <algorithm> // For std::remove_if

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
};