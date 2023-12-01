#pragma once // This is a common directive to prevent multiple inclusions of the header file

#include "utils.h"
#include <iostream>
#include <vector>
#include <queue>
#include <stack>

#define NONCORE_CLUSTER_FLAG -3

using namespace std;

struct MTNode
{
    MTNode() : epsilon_threshold(0.0), cluster_id(-1) {}
    MTNode(double epsilon, int c_id) : epsilon_threshold(epsilon), cluster_id(c_id) {}
    MTNode(double epsilon, const vector<int> &verts, int c_id) : epsilon_threshold(epsilon), vertices(verts), cluster_id(c_id) {}

    double epsilon_threshold;
    int cluster_id;
    vector<int> vertices; // vertex_id = vertex + query_type_offset
    vector<int> children;
};

class MuTree
{
public:
    vector<MTNode> MTNodes;
    int root_index;
    vector<float> CoreThresholdMu;
    // vector<int> leafIndices; // Stores indices of leaf nodes in MTNodes

    MuTree();

    void AddChildNode(int parent_index, double epsilon_threshold, const vector<int> &vertices, int c_id);
    void setMTNodes(int n);
    void FindandUpdateNode(int c_id, double epsilon, const vector<int> &vertices);
    int count_all_size();

    void save_tree(string save_path, int mu, int type_j);
    void load_tree(string load_path, int mu, int type_j);
    int query_tree(double e_threshold, vector<vector<int>> &res, vector<int> &core_id);

private:
    int FindNodeNumber(int c_id);

    void PrintFindAllNodeNumber(int c_id);
    void PrintDFSPaths(int start_index, vector<int> &current_path);
    void PrintDFSNode(int start_index);
};

// vector<MuTree> GetLeafNodes(const MuTree &root);
