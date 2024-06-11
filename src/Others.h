#pragma once

#include "HinGraph.h"
#include <vector>
#include <iostream>
#include <algorithm>

using namespace std;

class Others
{
public:
    Others(/* args */);
    ~Others();
    vector<bool> outlier;
    vector<bool> hub;
    vector<pair<int, int>> hub_cid;
    
    void compute_hub_outlier_index(const HinGraph &graph, const vector<int> &c_member_i, const vector<int> &community_number);
    void compute_hub_outlier_online(const HinGraph &graph, const vector<int> &c_member_i, const vector<int> &community_number);
    void output_CD_result(const HinGraph &graph, string output_path);

};
