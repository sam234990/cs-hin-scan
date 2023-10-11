#pragma once // This is a common directive to prevent multiple inclusions of the header file

#include "utils.h"
#include <vector>

using namespace std;

struct Query_nei_similarity
{
    Query_nei_similarity(int nei_id, float id_sim) : neighbor_id(nei_id), similarity(id_sim) {}

    int neighbor_id;
    float similarity;
};

struct OrderNeighbor
{
    vector<vector<Query_nei_similarity>> neighbor_;
};

struct ONdecouple
{
    vector<vector<int>> neighbor_id;
    vector<vector<float>> neighbor_sim;
};

bool compareBySimilarity(const Query_nei_similarity &a, const Query_nei_similarity &b);

class ONIndex
{
public:
    vector<OrderNeighbor> ONList;
    vector<ONdecouple> ONdelist;
    int all_node_size, type_size;

    ONIndex();
    void Initial_ON(int num_query_type, int type_num);

    void save_ON(string save_path, vector<int> type_j_vec);
    void load_ON(string load_path, vector<int> type_j_vec);
    // void query_ON(vector<int> type_order_vec);
};
