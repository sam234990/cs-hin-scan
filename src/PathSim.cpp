#include "PathSim.h"
#include "HinGraph.h"

PathSim::PathSim(/* args */)
{
}

PathSim::~PathSim()
{
}

std::string trim(const std::string &str)
{
    std::string trimmedStr = str;
    trimmedStr.erase(trimmedStr.begin(), std::find_if_not(trimmedStr.begin(), trimmedStr.end(), ::isspace));
    trimmedStr.erase(std::find_if_not(trimmedStr.rbegin(), trimmedStr.rend(), ::isspace).base(), trimmedStr.end());
    return trimmedStr;
}

void PathSim::initial_metapaths(vector<string> mps)
{
    MetaPathVec.resize(mps.size());
    for (int i = 0; i < mps.size(); i++)
    {
        string trimmedStr = trim(mps[i]);
        std::istringstream iss(trimmedStr);
        std::string token;

        // Split the string by spaces
        while (std::getline(iss, token, ' '))
        {
            int value = std::stoi(token);

            if (MetaPathVec[i].vertex.size() == MetaPathVec[i].edge.size())
                MetaPathVec[i].vertex.push_back(value);
            else
                MetaPathVec[i].edge.push_back(value);
        }
        MetaPathVec[i].pathLen = MetaPathVec[i].vertex.size();
        if (MetaPathVec[i].vertex.size() != (MetaPathVec[i].edge.size() + 1))
        {
            cerr << "wrong format metapath: " << mps[i] << endl;
            exit(-1);
        }
    }
}

void PathSim::initial_query_vertex(int q_num)
{
    query_num_ = q_num;
    query_path_cnt.resize(q_num, vector<pathcnt>(MetaPathVec.size()));
}

void PathSim::search(const HinGraph &graph, int query_vertex_id)
{
    int query_i = query_vertex_id - graph.query_type_offset_;
    for (size_t i = 0; i < MetaPathVec.size(); i++)
    {
        MetaPath &cur_p = MetaPathVec[i];
        queue<pair<int, int>> bfs_path;
        bfs_path.push(make_pair(query_vertex_id, 0));
        while (!bfs_path.empty())
        {
            int curVertex_id = bfs_path.front().first, cur_step = bfs_path.front().second;
            bfs_path.pop();
            if ((cur_step + 1) == cur_p.pathLen)
            {
                query_path_cnt[query_i][i].ins_path_cnt[curVertex_id - graph.query_type_offset_]++;
                continue;
            }

            // search its neighbor
            int nei_edge_start = graph.vertex_offset_[curVertex_id];
            int nei_end = ((curVertex_id + 1) == graph.n) ? graph.m : graph.vertex_offset_[curVertex_id + 1];
            for (int j = nei_edge_start; j < nei_end; j++)
            {
                int nei_type = graph.edges_[j].v_type, nei_id = graph.edges_[j].v_id, edge_type = graph.edges_[j].edge_type;
                if (nei_type != cur_p.vertex[cur_step + 1] || edge_type != cur_p.edge[cur_step])
                    continue; // cur_edge not follow the metapath
                bfs_path.push(make_pair(nei_id, cur_step + 1));
            }
        }
        query_path_cnt[query_i][i].finish_cnt = true;
    }
}

void PathSim::compute_pathsim(int i, int j, vector<float> &sim_res)
{
    sim_res.resize(MetaPathVec.size());
    for (size_t meta_i = 0; meta_i < MetaPathVec.size(); meta_i++)
    {
        int self_i = query_path_cnt[i][meta_i].ins_path_cnt[i];
        int self_j = query_path_cnt[j][meta_i].ins_path_cnt[j];
        int ins_ij = query_path_cnt[i][meta_i].ins_path_cnt[j];
        sim_res[meta_i] = float(ins_ij) / (self_i * self_j);
    }
}
