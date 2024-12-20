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
    path_num = mps.size();
    MetaPathVec.resize(path_num);
    OnlyVertexType.resize(path_num);
    for (int i = 0; i < path_num; i++)
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

    for (size_t i = 0; i < MetaPathVec.size(); i++)
    {
        MetaPath &cur_p = MetaPathVec[i];
        OnlyVertexType[i] = false;
        if (cur_p.edge[0] == -1)
            OnlyVertexType[i] = true;
    }
}

void PathSim::initial_query_vertex(int q_num)
{
    query_num_ = q_num;
    query_path_cnt.resize(q_num);
    path_search_finish_.resize(q_num);
    p_induce.resize(q_num);
    p_g.resize(q_num);
    for (int i = 0; i < q_num; i++)
    {
        path_search_finish_[i] = false;
        p_induce[i] = false;
        query_path_cnt[i] = vector<pathcnt>(path_num);
        p_g[i] = vector<int>();
    }
}

void PathSim::search(const HinGraph &graph, int query_i)
{
    if (path_search_finish_[query_i] == true)
        return;

    int query_vertex_id = query_i + graph.query_type_offset_;
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
                // if (nei_type != cur_p.vertex[cur_step + 1] || edge_type != cur_p.edge[cur_step])
                //     continue; // cur_edge not follow the metapath
                if (nei_type != cur_p.vertex[cur_step + 1])
                    continue;
                if ((cur_p.edge[cur_step] != -1) && (edge_type != cur_p.edge[cur_step]))
                    continue;
                bfs_path.push(make_pair(nei_id, cur_step + 1));
            }
        }
    }
    path_search_finish_[query_i] = true;
}

double PathSim::compute_avg_pathsim(const HinGraph &graph, int i, int j)
{
    if (i == j)
        return 1.0;
    if (path_search_finish_[i] == false)
        search(graph, i);
    if (path_search_finish_[j] == false)
        search(graph, j);
    double avg = 0.0;
    int num = 0;
    for (size_t meta_i = 0; meta_i < MetaPathVec.size(); meta_i++)
    {
        int self_i = query_path_cnt[i][meta_i].ins_path_cnt[i];
        int self_j = query_path_cnt[j][meta_i].ins_path_cnt[j];
        int ins_ij = query_path_cnt[i][meta_i].ins_path_cnt[j];
        double sim = static_cast<double>(2 * ins_ij) / (self_i + self_j);
        if ((self_i == 0) && (self_j == 0))
        {
            avg += 1.0;
            num++;
            continue;
        }
        avg += sim;
        num++;
    }
    return avg / num;
}

bool PathSim::judge_pathsim(int i, int j, const vector<float> sim_threshold)
{
    if (i == j)
        return true;
    for (size_t meta_i = 0; meta_i < MetaPathVec.size(); meta_i++)
    {
        int self_i = query_path_cnt[i][meta_i].ins_path_cnt[i];
        int self_j = query_path_cnt[j][meta_i].ins_path_cnt[j];
        int ins_ij = query_path_cnt[i][meta_i].ins_path_cnt[j];
        float sim = float(2 * ins_ij) / (self_i + self_j);
        if ((self_i == 0) && (self_j == 0))
            continue;
        if (sim < sim_threshold[meta_i])
            return false;
    }
    return true;
}

void PathSim::generate_cand_nei(const HinGraph &graph, int query_i, vector<int> &cand_nei)
{
    int num_cand = query_path_cnt[query_i][0].ins_path_cnt.size();
    int gen_cand_type = 0;
    for (int i = 0; i < path_num; i++)
    {
        if (num_cand > query_path_cnt[query_i][i].ins_path_cnt.size())
        {
            num_cand = query_path_cnt[query_i][i].ins_path_cnt.size();
            gen_cand_type = i;
        }
    }
    cand_nei.resize(num_cand);
    const auto &query_i_cnt = query_path_cnt[query_i][gen_cand_type].ins_path_cnt;
    int i = 0;
    for (auto it = query_i_cnt.begin(); it != query_i_cnt.end(); ++it, i++)
    {
        int nei_i = it->first;
        cand_nei[i] = nei_i + graph.query_type_offset_;
        if (path_search_finish_[nei_i] == false)
            search(graph, nei_i);
    }
    sort(cand_nei.begin(), cand_nei.end());
}

void PathSim::trans_homo_graph(const HinGraph &graph, string meta_path, string save_path)
{
    path_num = 1;
    MetaPath inst_path;
    string trimmedStr = trim(meta_path);
    std::istringstream iss(trimmedStr);
    std::string token;
    // Split the string by spaces
    while (std::getline(iss, token, ' '))
    {
        int value = std::stoi(token);

        if (inst_path.vertex.size() == inst_path.edge.size())
            inst_path.vertex.push_back(value);
        else
            inst_path.edge.push_back(value);
    }
    inst_path.pathLen = inst_path.vertex.size();
    if (inst_path.vertex.size() != (inst_path.edge.size() + 1))
    {
        cerr << "wrong format metapath: " << meta_path << endl;
        exit(-1);
    }

    cout << "finish load meta path" << endl;
    homo_graph.resize(graph.num_query_type_);
    homo_degree.resize(graph.num_query_type_);
    int m = 0;
    omp_set_num_threads(16);
#pragma omp parallel for reduction(+ : m) // 使用 reduction 以保证 m 的线程安全累加
    for (int i = 0; i < graph.num_query_type_; i++)
    {
        // homo_graph[i].reserve(1000);
        vector<int> i_neighbor;
        i_neighbor.reserve(1000);
        MetaPath &cur_p = inst_path;
        queue<pair<int, int>> bfs_path;
        bfs_path.push(make_pair(i + graph.query_type_offset_, 0));
        while (!bfs_path.empty())
        {
            int cur_id = bfs_path.front().first, cur_step = bfs_path.front().second;
            bfs_path.pop();
            if ((cur_step + 1) == cur_p.pathLen)
            {
                i_neighbor.push_back(cur_id - graph.query_type_offset_);
                continue;
            }

            // search its neighbor
            int nei_edge_start = graph.vertex_offset_[cur_id];
            int nei_end = ((cur_id + 1) == graph.n) ? graph.m : graph.vertex_offset_[cur_id + 1];
            for (int j = nei_edge_start; j < nei_end; j++)
            {
                int nei_type = graph.edges_[j].v_type, nei_id = graph.edges_[j].v_id, edge_type = graph.edges_[j].edge_type;
                // if (nei_type != cur_p.vertex[cur_step + 1] || edge_type != cur_p.edge[cur_step])
                //     continue; // cur_edge not follow the metapath
                if (nei_type != cur_p.vertex[cur_step + 1])
                    continue;
                if ((cur_p.edge[cur_step] != -1) && (edge_type != cur_p.edge[cur_step]))
                    continue;
                bfs_path.push(make_pair(nei_id, cur_step + 1));
            }
        }
        sort(i_neighbor.begin(), i_neighbor.end());
        i_neighbor.erase(unique(i_neighbor.begin(), i_neighbor.end()), i_neighbor.end());
        i_neighbor.erase(remove(i_neighbor.begin(), i_neighbor.end(), i), i_neighbor.end());
        int i_degree = i_neighbor.size();
        if (i_degree > graph.num_query_type_)
            cout << "error process this vertex, id: " << i + graph.query_type_offset_ << endl;
        m += i_degree;
        homo_degree[i] = i_degree;
        homo_graph[i] = move(i_neighbor);
        if (i % (graph.num_query_type_ / 10) == 0)
        {
#pragma omp critical
            cout << i << endl;
        }
    }

    cout << "finish trans homo graph, start save" << endl;
    string adj_file_path = save_path + "/adj.txt";
    ofstream adj_file = open_file_ofstream(adj_file_path);
    for (int i = 0; i < graph.num_query_type_; i++)
    {
        for (int id : homo_graph[i])
            adj_file << id << " ";
        adj_file << endl;
    }
    adj_file.close();

    string degree_file_path = save_path + "/degree.txt";
    ofstream degree_file = open_file_ofstream(degree_file_path);
    degree_file << "4\n"
                << graph.num_query_type_ << endl;
    degree_file << m << endl;
    for (int i = 0; i < graph.num_query_type_; i++)
    {
        degree_file << homo_degree[i] << endl;
    }
    degree_file.close();
    cout << "finish save" << endl;
}

vector<int> PathSim::p_induced_graph(const HinGraph &graph, int i)
{
    if (path_search_finish_[i] == false)
        search(graph, i);
    if (p_induce[i] == true)
        return p_g[i];

    vector<int> p_neighbor;
    int reserve_size = 0;
    for (size_t meta_i = 0; meta_i < MetaPathVec.size(); meta_i++)
        reserve_size += query_path_cnt[i][meta_i].ins_path_cnt.size();
    p_neighbor.reserve(reserve_size);
    for (size_t meta_i = 0; meta_i < MetaPathVec.size(); meta_i++)
    {
        for (const auto &pair : query_path_cnt[i][meta_i].ins_path_cnt)
            p_neighbor.push_back(pair.first);
    }
    sort(p_neighbor.begin(), p_neighbor.end());
    auto last = std::unique(p_neighbor.begin(), p_neighbor.end());
    p_neighbor.erase(last, p_neighbor.end());
    p_g[i] = move(p_neighbor);
    p_induce[i] = true;

    return p_g[i];
}

void PathSim::initial_hetesim(const HinGraph &graph, int vertex_num, int query_num)
{

    // HeteSim initial
    hete_path_cnt.resize(vertex_num);
    vertex_nei_type_split.resize(vertex_num);
    for (int i = 0; i < vertex_num; i++)
    {
        // HeteSim initial
        hete_path_cnt[i] = vector<path_hetesim>(path_num);
    }

    // get related vertex type
    set<int> vertex_type;
    for (size_t i = 0; i < MetaPathVec.size(); i++)
    {
        MetaPath &cur_p = MetaPathVec[i];
        for (int j = 0; j < cur_p.pathLen; j++)
            vertex_type.insert(cur_p.vertex[j]);
    }
    for (int type : vertex_type)
    {
        int type_start = graph.vertex_start_map_[type];
        int type_end = (type + 1) == graph.n_types ? graph.n : graph.vertex_start_map_[type + 1];
        for (int vertex_id = type_start; vertex_id < type_end; vertex_id++)
        {
            vertex_nei_type_split[vertex_id] = nei_list();
            int nei_start = graph.vertex_offset_[vertex_id];
            int nei_end = ((vertex_id + 1) == graph.n) ? graph.m : graph.vertex_offset_[vertex_id + 1];
            int nei_size = nei_end - nei_start;
            for (int j = nei_start; j < nei_end; j++)
            {
                int nei_type = graph.edges_[j].v_type, nei_id = graph.edges_[j].v_id;
                int edge_type = graph.edges_[j].edge_type;
                // add vertex nei list
                if (vertex_nei_type_split[vertex_id].vertex_type_nei[nei_type].empty())
                    vertex_nei_type_split[vertex_id].vertex_type_nei[nei_type].reserve(nei_size);
                vertex_nei_type_split[vertex_id].vertex_type_nei[nei_type].push_back(nei_id);

                // add edge type nei list
                if (vertex_nei_type_split[vertex_id].edge_type_nei[edge_type].empty())
                    vertex_nei_type_split[vertex_id].edge_type_nei[edge_type].reserve(nei_size);
                vertex_nei_type_split[vertex_id].edge_type_nei[edge_type].push_back(nei_id);
            }
        }
    }
}

bool PathSim::judge_hetesim(int id1, int id2, const vector<float> sim_threshold)
{
    if (id1 == id2)
        return true;
    for (size_t meta_i = 0; meta_i < MetaPathVec.size(); meta_i++)
    {
        double sim = compute_hetesim(id1, id2, meta_i, 0);
        if (sim < sim_threshold[meta_i])
            return false;
    }
    return true;
}

double PathSim::compute_hetesim(int id1, int id2, int mp_id, int step)
{
    if (id1 == id2)
        return 1.0;
    if (step == MetaPathVec[mp_id].pathLen / 2) // only should compute half of the path
        return 0.0;

    if (hete_path_cnt[id1][mp_id].hetesim.find(id2) != hete_path_cnt[id1][mp_id].hetesim.end())
        return hete_path_cnt[id1][mp_id].hetesim[id2];
    if (hete_path_cnt[id2][mp_id].hetesim.find(id1) != hete_path_cnt[id2][mp_id].hetesim.end())
        return hete_path_cnt[id2][mp_id].hetesim[id1];

    // compute hetesim
    int mp_len = MetaPathVec[mp_id].pathLen;
    std::vector<int> *id1_nei_list = nullptr;
    std::vector<int> *id2_nei_list = nullptr;
    if (OnlyVertexType[mp_id] == true)
    {
        int next_type_1 = MetaPathVec[mp_id].vertex[step + 1];
        int next_type_2 = MetaPathVec[mp_id].vertex[mp_len - step - 2];
        id1_nei_list = &vertex_nei_type_split[id1].vertex_type_nei[next_type_1];
        id2_nei_list = &vertex_nei_type_split[id2].vertex_type_nei[next_type_2];
    }
    else
    {
        int next_type_1 = MetaPathVec[mp_id].edge[step];
        int next_type_2 = MetaPathVec[mp_id].edge[mp_len - step - 2];
        id1_nei_list = &vertex_nei_type_split[id1].edge_type_nei[next_type_1];
        id2_nei_list = &vertex_nei_type_split[id2].edge_type_nei[next_type_2];
    }

    int id1_size = id1_nei_list->size();
    int id2_size = id2_nei_list->size();

    if (id1_size == 0 || id2_size == 0)
    {
        hete_path_cnt[id1][mp_id].hetesim[id2] = 0.0;
        hete_path_cnt[id2][mp_id].hetesim[id1] = 0.0;
        return 0.0;
    }

    double sum_sim = 0.0;
    for (auto id1_nei : *id1_nei_list)
        for (auto id2_nei : *id2_nei_list)
            sum_sim += compute_hetesim(id1_nei, id2_nei, mp_id, step + 1);
    
    double denominator = sqrt(id1_size * id2_size);
    double sim = sum_sim / denominator;
    hete_path_cnt[id1][mp_id].hetesim[id2] = sim;
    hete_path_cnt[id2][mp_id].hetesim[id1] = sim;
    return sim;
}
