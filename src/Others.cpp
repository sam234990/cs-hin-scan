#include "Others.h"
#include "HinGraph.h"

#define OUTLIER_INIT -1

Others::Others(/* args */)
{
}

Others::~Others()
{
}

void Others::compute_hub_outlier_index(const HinGraph &graph, const vector<int> &c_member_i, const vector<int> &community_number)
{
    vector<bool> visit_mem(graph.num_query_type_, false);
    vector<int> outlier_cid(graph.num_query_type_, OUTLIER_INIT);
    hub_cid.reserve(c_member_i.size());

    hub.resize(graph.num_query_type_);
    outlier.resize(graph.num_query_type_);

    for (const int mem_i : c_member_i)
    {
        if (graph.community_member[mem_i] == false)
            continue;
        if (visit_mem[mem_i] == true)
            continue;
        visit_mem[mem_i] = true;
        for (const auto &i_nei : graph.h_sim[mem_i])
        {
            if (graph.community_member[i_nei.neighbor_i] == true)
                continue; // vertex in SSC
            if (outlier_cid[i_nei.neighbor_i] == OUTLIER_INIT)
            { // outlier vertex
                outlier[i_nei.neighbor_i] = true;
                outlier_cid[i_nei.neighbor_i] = community_number[mem_i];
            }
            else
            { // outlier vertex
                if (outlier_cid[i_nei.neighbor_i] != community_number[mem_i])
                {
                    hub[i_nei.neighbor_i] = true;
                    hub_cid.push_back(make_pair(i_nei.neighbor_i, community_number[mem_i]));
                }
            }
        }
    }
}

void Others::compute_hub_outlier_online(const HinGraph &graph, const vector<int> &c_member_i, const vector<int> &community_number)
{
    vector<bool> visit_mem(graph.num_query_type_, false);
    vector<int> outlier_cid(graph.num_query_type_, OUTLIER_INIT);

    hub_cid.reserve(c_member_i.size());

    hub.resize(graph.num_query_type_);
    outlier.resize(graph.num_query_type_);

    for (const int mem_i : c_member_i)
    {
        if (graph.community_member[mem_i] == false)
            continue;
        if (visit_mem[mem_i] == true)
            continue;
        visit_mem[mem_i] = true;
        for (const auto &nei_id : graph.cand_sn_list[mem_i])
        {
            int nei_i = nei_id - graph.query_type_offset_;
            if (graph.community_member[nei_i] == true)
                continue; // vertex in SSC
            if (outlier_cid[nei_i] == OUTLIER_INIT)
            { // outlier vertex
                outlier[nei_i] = true;
                outlier_cid[nei_i] = community_number[mem_i];
            }
            else
            { // hub vertex
                if (outlier_cid[nei_i] != community_number[mem_i])
                {
                    hub[nei_i] = true;
                    hub_cid.push_back(make_pair(nei_i, community_number[mem_i]));
                }
            }
        }
    }
}

void Others::output_CD_result(const HinGraph &graph, string output_path)
{
    string output_file_path = output_path + "/cd" + to_string(graph.p_query_type) + "-mu" + to_string(graph.p_mu);
    output_file_path += "-" + graph.query_file_name + ".txt";
    ofstream output_file = open_file_ofstream(output_file_path);
    cout << "output file name and path: " << output_file_path << endl;
    cout << "vertex_id -- cluster_id\n";

    for (const auto &res : graph.all_res_com)
    {
        int cid = -1;
        for (const auto ci : res)
        {
            if (graph.cand_core_[ci] == true)
            {
                cid = ci;
                break;
            }
        }
        if (cid == -1)
            continue;

        for (const auto ci : res)
        {
            int vertex_id = ci + graph.query_type_offset_;
            output_file << vertex_id << " " << cid << endl;
        }
    }

    output_file.close();
}