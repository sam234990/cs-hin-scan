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
    bool finish_cnt;
};

class HinGraph; 

class PathSim
{
public:
    std::vector<MetaPath> MetaPathVec;
    int query_num_;
    vector<vector<pathcnt>> query_path_cnt;

    PathSim(/* args */);
    ~PathSim();
    void initial_metapaths(vector<string> mps);
    void initial_query_vertex(int q_num);
    void search(const HinGraph& graph, int query_vertex_id);
    void compute_pathsim(int i, int j, vector<float> &sim_res);
};

