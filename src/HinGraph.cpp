#include "HinGraph.h"

int initial_vector_size = 32;

void computeIntersection_set(vector<int> &vec1, const vector<int> &vec2)
{
    vector<int> intersection;
    int i = 0; // pointer for vec1
    int j = 0; // pointer for vec2

    while (i < vec1.size() && j < vec2.size())
    {
        if (vec1[i] == vec2[j])
        {
            intersection.push_back(vec1[i]);
            i++;
            j++;
        }
        else if (vec1[i] < vec2[j])
        {
            i++;
        }
        else
        {
            j++;
        }
    }
    vec1 = move(intersection);
    return;
}

void computeIntersection_QN(vector<Query_nei_similarity> &vec1, const vector<Query_nei_similarity> &vec2)
{
    vector<Query_nei_similarity> intersection;
    intersection.reserve(vec1.size());
    int i = 0; // 指向 vec1 的指针
    int j = 0; // 指向 vec2 的指针

    while (i < vec1.size() && j < vec2.size())
    {
        if (vec1[i].neighbor_id == vec2[j].neighbor_id)
        {
            intersection.push_back(vec1[i]);
            i++;
            j++;
        }
        else if (vec1[i].neighbor_id < vec2[j].neighbor_id)
        {
            i++;
        }
        else
        {
            j++;
        }
    }
    vec1 = move(intersection);
    return;
}

void computeIntersection_QN2(vector<Query_nei_similarity> &vec1, const vector<int> &vec2)
{
    vector<Query_nei_similarity> intersection;
    intersection.reserve(vec2.size());
    int i = 0; // 指向 vec1 的指针
    int j = 0; // 指向 vec2 的指针

    while (i < vec1.size() && j < vec2.size())
    {
        if (vec1[i].neighbor_id == vec2[j])
        {
            intersection.push_back(vec1[i]);
            i++;
            j++;
        }
        else if (vec1[i].neighbor_id < vec2[j])
        {
            i++;
        }
        else
        {
            j++;
        }
    }
    vec1 = move(intersection);
    return;
}

vector<int> randomSampling(int n, int sampleSize)
{
    vector<int> samples;
    int interval = n / sampleSize; // internal
    for (int i = 1; i < n; i += interval)
    {
        samples.push_back(i);
    }
    return samples;
}

HinGraph::HinGraph(string data_dir)
{
    size_t lastSlashPos = data_dir.find_last_of('/');
    if (lastSlashPos != std::string::npos)
    {
        data_file_name = data_dir.substr(lastSlashPos + 1);
    }
    else
    {
        data_file_name = data_dir;
    }

    data_dir_ = data_dir + "/pro_graph";
    data_index_dir_ = data_dir + "/index";
}

HinGraph::~HinGraph()
{
    if (edges_ != NULL)
        delete[] edges_;
    if (vertex_offset_ != NULL)
        delete[] vertex_offset_;
}

void HinGraph::load_graph()
{
    cout << "start load graph" << endl;
    string graph_info_path = data_dir_ + "/graph_info.txt";
    string line;
    ifstream graph_info = open_file_fstream(graph_info_path);
    getline(graph_info, line);
    istringstream iss(line);
    iss >> m >> n >> n_types;
    hin_schema_adjacencyMatrix.resize(n_types, ::vector<int>(n_types, 0));
    getline(graph_info, line);
    istringstream iss2(line);
    for (int i = 0; i < n_types; i++)
    {
        int start_num;
        iss2 >> start_num;
        vertex_start_map_.push_back(start_num);
    }
    for (int i = 0; i < n_types; i++)
    {
        getline(graph_info, line);
        istringstream iss3(line);

        for (int j = 0; j < n_types; j++)
        {
            iss3 >> hin_schema_adjacencyMatrix[i][j];
        }
    }
    for (int i = 0; i < n_types; i++)
    {
        getline(graph_info, line);
        istringstream iss4(line);
        for (int j = 0; j < n_types; j++)
        {
            iss4 >> hin_schema_edge_cnt[i][j];
        }
    }
    graph_info.close();

    cout << m << " " << n << " " << n_types << endl;
    vertex_offset_ = new int[n];
    edges_ = new Vertex_type[m];

    string graph_type_path = data_dir_ + "/graph_type.txt";
    ifstream graph_type = open_file_fstream(graph_type_path);
    for (int i = 0; i < m; i++)
    {
        getline(graph_type, line);
        istringstream iss_graph_type(line);
        iss_graph_type >> edges_[i].v_id >> edges_[i].v_type;
    }
    graph_type.close();

    string vertex_offset_path = data_dir_ + "/graph_offset.txt";
    ifstream vertex_offset = open_file_fstream(vertex_offset_path);
    for (int i = 0; i < n; i++)
    {
        getline(vertex_offset, line);
        istringstream iss_vertex_offset(line);
        iss_vertex_offset >> vertex_offset_[i];
    }
    cout << "finish read graph" << endl;
    vertex_offset.close();
    return;
}

void HinGraph::cs_hin_scan(string query_file, string mode)
{
    load_query_file(query_file);
    initialize_query_();
    reinitialize_query_();

    if (mode == "-qidx" || mode == "-qidx1" || mode == "-qidxON")
    {
        if (mode == "-qidx1")
        {
            mode_query = 11;
        }
        if (mode == "-qidxON")
        {
            mode_query = 111;
        }
        return;
    }
    else
    {
        if (mode == "-q")
        {
            baseline_query_();
        }

        if (mode == "-q1")
        {
            mode_query = 1;
        }
        else if (mode == "-q2")
        {
            mode_query = 2;
        }
        else if (mode == "-q3")
        {
            mode_query = 3;
        }
        // improved_query_();
    }

    // count_core_vertices();
}

void HinGraph::load_query_file(string query_file_path)
{
    size_t lastSlashPos = query_file_path.find_last_of('/');
    if (lastSlashPos != string::npos)
    {
        string fileName = query_file_path.substr(lastSlashPos + 1);
        size_t lastDotPos = fileName.find_last_of('.');
        if (lastDotPos != string::npos)
        {
            query_file_name = fileName.substr(0, lastDotPos);
        }
    }
    else
    {
        query_file_name = query_file_path;
    }
    cout << "start read query paramets" << endl;
    unit_epsilon = false;
    string line;
    ifstream query_file = open_file_fstream(query_file_path);
    int lineCount = 0;
    random_query = false;
    while (getline(query_file, line))
    {
        istringstream iss(line);
        if (lineCount == 0)
        { // 1. mu
            if (!(iss >> p_mu))
            {
                cerr << "Failed to read mu from file." << endl;
                exit(1);
            }
            if (p_mu <= 0)
            {
                cerr << "Parameter mu setting Error." << p_mu << endl;
                exit(1);
            }
        }
        else if (lineCount == 1)
        { // 2. distance
            if (!(iss >> p_d))
            {
                cerr << "Failed to read distance from file." << endl;
                exit(1);
            }
            if (p_d < 0 || p_d >= n_types)
            {
                cerr << "Parameter distance setting Error." << p_d << endl;
                exit(1);
            }
        }
        else if (lineCount == 2)
        { // 3. query_type
            if (!(iss >> p_query_type >> query_node_num))
            {
                cerr << "Failed to read query_type from file." << endl;
                exit(1);
            }
            if (p_query_type < 0 || p_query_type >= n_types)
            {
                cerr << "Parameter query type setting Error." << p_query_type << endl;
                exit(1);
            }
            query_node_list.reserve(query_node_num);
        }
        else if (lineCount == 3)
        { // 4. query node list
            for (int i = 0; i < query_node_num; i++)
            {
                iss >> query_node_list[i];
                if (i == 0 && query_node_list[i] == -1)
                {
                    random_query = true;
                    break;
                }
            }
        }
        else
        { // 5. epsilon map
            int key;
            double value;
            if (iss >> key >> value)
            {
                if (key == -1)
                {
                    unit_epsilon = true;
                    unit_epsilon_value = value;
                    cout << "Use the unit epsilon vector" << endl;
                    break;
                }
                else
                    type_epsilon[key] = value;
            }
            else
            {
                cerr << "Failed to read map entry from file." << endl;
                exit(1);
            }
            if (key < 0 || key >= n_types)
            {
                cerr << "Parameter other_type setting Error." << key << endl;
                exit(1);
            }
            if (value < 0 || value > 1)
            {
                cerr << "Parameter epsilon setting Error." << value << endl;
                exit(1);
            }
        }
        lineCount++;
    }
    query_file.close();

    // 输出读取结果
    cout << "mu: " << p_mu << endl;
    cout << "distance: " << p_d << endl;
    cout << "query_type: " << p_query_type << " query_node_num: " << query_node_num << endl;
    cout << "type_epsilon: ";
    if (unit_epsilon)
    {
        cout << unit_epsilon_value << endl;
    }
    else
    {
        for (const auto &pair : type_epsilon)
        {
            cout << pair.first << " -> " << pair.second << " -- ";
        }
        cout << endl;
    }
    return;
}

void HinGraph::initialize_query_()
{
    vector<bool> visited(n_types, false);
    queue<int> q;

    // get schema distance
    q.push(p_query_type);
    visited[p_query_type] = true;
    distance_[p_query_type] = 0;
    while (!q.empty())
    {
        int currVertex = q.front();
        q.pop();
        int currDistance = distance_[currVertex];
        for (int neighbor = 0; neighbor < n_types; ++neighbor)
        {
            if (hin_schema_adjacencyMatrix[currVertex][neighbor] == 1 && !visited[neighbor])
            {
                q.push(neighbor);
                visited[neighbor] = true;
                distance_[neighbor] = currDistance + 1;
            }
        }
    }

    if (unit_epsilon)
    {
        for (const auto &entry : distance_)
        {
            int vertex_type = entry.first, distance = entry.second;
            if (distance > p_d)
                continue;
            if (p_d == 1 && distance == 0)
                continue;
            type_epsilon[vertex_type] = unit_epsilon_value;
            cout << vertex_type << " -> " << type_epsilon[vertex_type] << " -- ";
            if (distance >= 1)
            {
                int type_start = vertex_start_map_[vertex_type];
                int type_end = (vertex_type + 1) == n_types ? n : vertex_start_map_[vertex_type + 1];
                for (int j = type_start; j < type_end; j++)
                {
                    t_hop_visited[j] = vector<int>();
                    t_hop_visited.reserve(initial_vector_size);
                }
            }
        }
    }
    else
    {
        for (const auto &entry : distance_)
        {
            int vertex = entry.first;
            int distance = entry.second;
            if (vertex != p_query_type && type_epsilon.count(vertex) == 0)
            { // this type neighbor is unconsidered.
                distance_[vertex] = -1;
                distance = -1;
            }
            // cout << "Shortest distance from vertex " << p_query_type << " to vertex " << vertex << ": " << distance_[vertex] << endl;
            if (distance >= 1)
            {
                int type_start = vertex_start_map_[vertex];
                int type_end = (vertex + 1) == n_types ? n : vertex_start_map_[vertex + 1];
                for (int j = type_start; j < type_end; j++)
                { // create a visited set for each two neighbor
                    t_hop_visited[j] = vector<int>();
                    t_hop_visited.reserve(initial_vector_size);
                }
                // cout << "create a visited unordered_set for each two neighbor :";
                // cout << type_start << "-" << type_end << endl;
            }
        }
    }
    int min_distance = 1000;
    for (auto &pair : type_epsilon)
    {
        if (pair.first == p_query_type)
            continue;
        if (pair.second == 0.0)
            continue;
        if (distance_[pair.first] < min_distance)
        {
            cand_gen_type = pair.first;
            min_distance = distance_[pair.first];
        }
    }

    // initial query neighbor
    query_type_offset_ = vertex_start_map_[p_query_type];
    int end = ((p_query_type + 1) == n_types) ? n : vertex_start_map_[p_query_type + 1];
    num_query_type_ = end - query_type_offset_;
    // cout << num_query_type_ << endl;

    if (random_query)
    {
        query_node_list = move(randomSampling(num_query_type_, query_node_num));
    }
    // initial type neighbor
    dn_adj_List.resize(num_query_type_);
    qn_adj_List.resize(num_query_type_);
    cand_sn_list.resize(num_query_type_);
    non_sn_list.resize(num_query_type_);
    visited_qtv_.resize(num_query_type_);
    ks_visit.resize(num_query_type_);
    is_in_community.resize(num_query_type_);
    similar_degree.resize(num_query_type_);
    effective_degree.resize(num_query_type_);
    cand_core_.resize(num_query_type_);
    cout << "finish initial" << endl;
    return;
}

void HinGraph::reinitialize_query_()
{
    for (int i = 0; i < num_query_type_; i++)
    {
        for (const auto &pair : type_epsilon)
        {
            dn_adj_List[i].d_neighbor_[pair.first] = vector<int>();
            dn_adj_List[i].d_neighbor_[pair.first].reserve(initial_vector_size);
        }
        qn_adj_List[i] = vector<int>();
        qn_adj_List[i].reserve(initial_vector_size);
        cand_sn_list[i] = vector<int>();
        cand_sn_list[i].reserve(initial_vector_size);
        non_sn_list[i] = vector<int>();
        non_sn_list[i].reserve(initial_vector_size);
        visited_qtv_[i] = false;
        ks_visit[i] = false;
        is_in_community[i] = false;
        similar_degree[i] = 0;
        effective_degree[i] = 0;
        cand_core_[i] = false;
    }
    // cout << "finish reinitialize variable" << endl;
    return;
}

void HinGraph::print_result()
{
    int num_community = 0, core_num = 0;
    for (int i = 0; i < num_query_type_; i++)
    {
        if (is_in_community[i])
            num_community++;
        if (cand_core_[i])
            core_num++;
    }
    if (similar_degree[query_i] >= p_mu)
        cout << query_i + query_type_offset_ << " in community contains number of vertex is : " << num_community << " Core number : " << core_num << endl;
}

void HinGraph::baseline_query_()
{
    cout << "start baseline online query " << endl;
    long all_time = 0;
    vector<long> time_cost(query_node_num, 0);
    for (int i = 0; i < query_node_num; i++)
    {
        query_i = query_node_list[i];
        reinitialize_query_();
        queue<int> cs_node_q;
        is_in_community[query_i] = true;

        Timer t1;
        t1.Start();
        search_k_strata(query_i);
        search_cand_sn(query_i);
        cs_check_cluster_core(query_i, cs_node_q);
        if (similar_degree[query_i] < p_mu)
        {
            cout << query_i << " Cannot search a community" << endl;
        }
        else
        {
            explore_community(cs_node_q);
        }

        // all_time += t1.StopAndPrint("finish baseline online query. Use time");
        long cost_time = t1.StopTime();
        all_time += cost_time;
        time_cost[i] = cost_time;
        print_result();
    }
    string str1 = "finish query " + to_string(query_node_num) + " times, use time:";
    Timer::PrintTime(str1, all_time);
    // cout << "finish query " << query_node_num << " times, use time: " << all_time << endl;
}

// void HinGraph::estimate_best_order(vector<int> ){}

void HinGraph::improved_query_()
{

    cout << "start improve online query " << endl;
    long all_time = 0;
    vector<long> time_cost(query_node_num, 0);
    for (int i = 0; i < query_node_num; i++)
    {
        query_i = query_node_list[i];
        reinitialize_query_();
        queue<int> cs_node_q;
        is_in_community[query_i] = true;

        Timer t1;
        t1.Start();
        search_k_strata(query_i);
        search_cand_sn(query_i);
        cs_check_cluster_core(query_i, cs_node_q);
        if (similar_degree[query_i] < p_mu)
        {
            cout << query_i << " Cannot search a community" << endl;
        }
        else
        {
            explore_community(cs_node_q);
        }

        long cost_time = t1.StopTime();
        all_time += cost_time;
        time_cost[i] = cost_time;
        print_result();
    }
    string str1 = "finish query " + to_string(query_node_num) + " times, use time:";
    Timer::PrintTime(str1, all_time);
}

void HinGraph::search_k_strata(int i)
{
    int query_vertex_id = i + query_type_offset_;
    // cout << "query_vertex_id: " << i << " \n";
    queue<pair<int, int>> bfs_path;
    unordered_set<int> v_vertex;
    bfs_path.push(make_pair(query_vertex_id, 0));
    v_vertex.insert(query_vertex_id);
    while (!bfs_path.empty())
    {
        int curVertex_id = bfs_path.front().first, cur_step = bfs_path.front().second;
        int cur_type = get_vertex_type(curVertex_id);
        // int cur_dis = distance_[cur_type];
        // cout << curVertex_id << " " << cur_type << " " << cur_dis << " -- ";
        bfs_path.pop();

        dn_adj_List[i].d_neighbor_[cur_type].push_back(curVertex_id);
        if (cur_step >= p_d)
        { // search finish
            continue;
        }
        // judge its neighbor
        int nei_edge_start = vertex_offset_[curVertex_id];
        int nei_end = ((curVertex_id + 1) == n) ? m : vertex_offset_[curVertex_id + 1];

        for (int j = nei_edge_start; j < nei_end; j++)
        {
            int nei_type = edges_[j].v_type, nei_id = edges_[j].v_id;
            if (v_vertex.find(nei_id) != v_vertex.end())
            { // this neighbor visited before
                continue;
            }
            int nei_dis = distance_[nei_type];
            if (nei_dis == -1)
            { // this type neighbor is unconsidered. continue
                continue;
            }

            v_vertex.insert(nei_id);
            bfs_path.push(make_pair(nei_id, cur_step + 1));
        }
    }
    for (const auto &pair : type_epsilon)
    {
        int type = pair.first;
        double type_epsilon = pair.second;
        if (type_epsilon == 0.0)
        {
            continue;
        }
        dn_adj_List[i].num_d_neighbor[type] = dn_adj_List[i].d_neighbor_[type].size();
        sort(dn_adj_List[i].d_neighbor_[type].begin(), dn_adj_List[i].d_neighbor_[type].end());
    }
    ks_visit[i] = true;
}

void HinGraph::search_cand_sn(int i)
{
    // 1. target type is considered and located in k-strata
    Vertex_neighbor &dn_i = dn_adj_List[i];
    // if (distance_[p_query_type] != -1 && dn_i.num_d_neighbor[p_query_type] != 0)
    // {
    //     cand_sn_list[i].assign(dn_i.d_neighbor_[p_query_type].begin(), dn_i.d_neighbor_[p_query_type].end());
    //     return;
    // }

    // 2. target type is not considered, use cand gen type k-strata search cand_sn
    int cand_gen_dis = distance_[cand_gen_type];
    const vector<int> &gen_type_dn = dn_i.d_neighbor_[cand_gen_type];
    unordered_set<int> v_vertex;
    for (auto v_i : gen_type_dn)
    {
        v_vertex.insert(v_i);
        queue<pair<int, int>> gen_q;
        gen_q.push(make_pair(v_i, 0));
        while (!gen_q.empty())
        {
            int vertex_i = gen_q.front().first, cur_step = gen_q.front().second;
            int cur_type = get_vertex_type(vertex_i);
            int cur_dis = distance_[cur_type];
            gen_q.pop();
            if (cur_step >= cand_gen_dis)
                continue;
            int nei_edge_start = vertex_offset_[vertex_i];
            int nei_end = ((vertex_i + 1) == n) ? m : vertex_offset_[vertex_i + 1];
            for (int j = nei_edge_start; j < nei_end; j++)
            {
                int nei_type = edges_[j].v_type, nei_id = edges_[j].v_id;
                if (nei_type == p_query_type)
                {
                    if (v_vertex.find(nei_id) != v_vertex.end()) // this neighbor visited before
                        continue;

                    v_vertex.insert(nei_id);
                    cand_sn_list[i].push_back(nei_id);
                }
                else
                {
                    if (distance_[nei_type] > cur_dis) // wrong direction
                        continue;
                    if (v_vertex.find(nei_id) != v_vertex.end()) // this neighbor visited before
                        continue;

                    v_vertex.insert(nei_id);
                    gen_q.push(make_pair(nei_id, cur_step + 1));
                }
            }
        }
    }
    sort(cand_sn_list[i].begin(), cand_sn_list[i].end());
}

void HinGraph::cs_check_cluster_core(int u, queue<int> &cs_queue)
{
    // 1. check core
    int vertex_u_id = u + query_type_offset_;
    int j = 0;
    unordered_set<int> v_vertex(qn_adj_List[u].begin(), qn_adj_List[u].end());
    for (auto not_i : non_sn_list[u])
        v_vertex.insert(not_i);
    if (similar_degree[u] < p_mu)
    {
        similar_degree[u] = qn_adj_List[u].size();
        effective_degree[u] = cand_sn_list[u].size();
        for (; j < cand_sn_list[u].size(); j++)
        {
            int v_id = cand_sn_list[u][j];
            int v = v_id - query_type_offset_;
            if (ks_visit[v] == false) // this vertex has not searched k-strata
                search_k_strata(v);
            if (v_vertex.find(v_id) != v_vertex.end()) // this sn has been computed before
                continue;

            bool sim_res = check_struc_sim(u, v);
            if (sim_res)
            {
                similar_degree[u]++;
                qn_adj_List[u].push_back(v_id);
            }
            else
                effective_degree[u]--;

            if (visited_qtv_[v] == false && v != u)
            {
                if (sim_res)
                {
                    similar_degree[v]++;
                    qn_adj_List[v].push_back(vertex_u_id);
                }
                else
                    non_sn_list[v].push_back(vertex_u_id);
            }
            if (effective_degree[u] < p_mu || similar_degree[u] >= p_mu)
                break;
        }
    }
    visited_qtv_[u] = true;

    // 2. cluster core
    if (similar_degree[u] < p_mu)
    {
        cand_core_[u] = false;
        return;
    }
    cand_core_[u] = true;
    for (auto visit_sn_id : qn_adj_List[u])
    {
        int sn_i = visit_sn_id - query_type_offset_;
        if (is_in_community[sn_i] == false)
        {
            is_in_community[sn_i] = true; // add into community
            cs_queue.push(sn_i);          // add to expand queue
        }
    }
    for (; j < cand_sn_list[u].size(); j++)
    {
        int v_id = cand_sn_list[u][j];
        int v = v_id - query_type_offset_;
        if (is_in_community[v] == true) // already in community
            continue;
        if (v_vertex.find(v_id) != v_vertex.end()) // this sn has been computed before
            continue;
        if (ks_visit[v] == false) // this vertex has not searched k-strata
            search_k_strata(v);

        bool sim_res = check_struc_sim(u, v);
        if (sim_res)
        {
            similar_degree[u]++;
            qn_adj_List[u].push_back(v_id);
            if (is_in_community[v] == false)
            {
                is_in_community[v] = true;
                cs_queue.push(v); // add to expand queue
            }
        }

        if (visited_qtv_[v] == false && v != u)
        {
            if (sim_res)
            {
                similar_degree[v]++;
                qn_adj_List[v].push_back(vertex_u_id);
            }
            else
                non_sn_list[v].push_back(vertex_u_id);
        }
    }
}

void HinGraph::explore_community(queue<int> &cs_queue)
{
    while (!cs_queue.empty())
    {
        int vertex_i = cs_queue.front();
        cs_queue.pop();
        search_cand_sn(vertex_i);
        cs_check_cluster_core(vertex_i, cs_queue);
    }
}

void computeIntersection_cnt(vector<int> &vec1, const vector<int> &vec2, vector<vector<int>> &qn_cnt,
                             const vector<int> &cnt2, int type_i)
{
    vector<int> intersection;
    vector<vector<int>> intersection_cnt;
    int i = 0; // 指向 vec1 的指针
    int j = 0; // 指向 vec2 的指针

    while (i < vec1.size() && j < vec2.size())
    {
        if (vec1[i] == vec2[j])
        {
            intersection.push_back(vec1[i]);
            qn_cnt[i].push_back(cnt2[j]);
            intersection_cnt.push_back(qn_cnt[i]);
            i++;
            j++;
        }
        else if (vec1[i] < vec2[j])
        {
            i++;
        }
        else
        {
            j++;
        }
    }
    vec1 = move(intersection);
    qn_cnt = move(intersection_cnt);
    // return intersection;
    return;
}

void HinGraph::check_neighbor(int i)
{
    const Vertex_neighbor &vertex = dn_adj_List[i];
    for (const auto &kvp : vertex.d_neighbor_)
    {
        int key = kvp.first;
        const vector<int> &values = kvp.second;

        set<int> uniqueValues(values.begin(), values.end());

        if (values.size() != uniqueValues.size())
        {
            cout << "Vector in element " << i << " with key " << key << " contains duplicate numbers. " << values.size() << "--" << uniqueValues.size() << endl;
        }
        else
        {
            cout << "Vector in element " << i << " with key " << key << " does not contain duplicate numbers." << values.size() << endl;
        }
    }
    return;
    const vector<int> &qn_adj = qn_adj_List[i];
    cout << i << " th Query Neighborhood size is : " << qn_adj.size() << " first 10 Neighbor : " << endl;
    int count = 0;
    bool contain_self = false;
    for (auto it = qn_adj.begin(); it != qn_adj.end() && count < 10; ++it, ++count)
    {
        if (count < 10)
            cout << *it << " ";
        if (*it == (i + query_type_offset_))
            contain_self = true;
    }
    cout << endl;
    if (contain_self == false)
    {
        cout << "this vertex QN don't contain itself" << endl;
        exit(-1);
    }
}

void HinGraph::check_two_hop_visited()
{
    for (const auto &kvp : t_hop_visited)
    {
        int vertex_id = kvp.first;
        int vertex_type = get_vertex_type(vertex_id);
        const vector<int> &visited = kvp.second;
        set<int> uniqueValues(visited.begin(), visited.end());

        if (visited.size() != uniqueValues.size())
        {
            cout << "Vertex id: " << vertex_id << " Vertex type: " << vertex_type << " contains duplicate visited vertices. " << visited.size() << "--" << uniqueValues.size() << endl;
        }
        // else
        // {
        //     cout << "Vertex id: " << vertex_id << " Vertex type: " << vertex_type << " does not contain duplicate numbers." << visited.size() << endl;
        // }
    }
}

bool HinGraph::check_struc_sim(int a, int b)
{
    if (a == b)
        return true;

    for (const auto &pair : type_epsilon)
    {
        int type = pair.first;
        double type_epsilon = pair.second;
        if (type_epsilon == 0.0)
            continue;

        if (judgeJacSim(dn_adj_List[a].d_neighbor_[type], dn_adj_List[b].d_neighbor_[type], type_epsilon) == false)
            return false;
    }
    return true;
}

bool HinGraph::judgeJacSim(const vector<int> &vec1, const vector<int> &vec2, double type_i_epsilon)
{
    int size_a = vec1.size(), size_b = vec2.size();
    if (size_a == 0 && size_b == 0)
    {
        return true;
    }
    else if ((size_a == 0 && size_b != 0) || (size_a != 0 && size_b == 0))
    {
        return false;
    }
    int intersection = 0;
    size_t i = 0;
    size_t j = 0;

    while (i < vec1.size() && j < vec2.size())
    {
        if (vec1[i] == vec2[j])
        {
            intersection++;
            i++;
            j++;
        }
        else if (vec1[i] < vec2[j])
        {
            i++;
        }
        else
        {
            j++;
        }
    }
    int unionSize = size_a + size_b - intersection;
    double sim = static_cast<double>(intersection) / static_cast<double>(unionSize);

    return sim >= type_i_epsilon;
}

void HinGraph::count_core_vertices()
{
    cout << "count all core vertices" << endl;
    int core_num = 0, print_num = 0;

    for (int u = 0; u < num_query_type_; u++)
    {
        int vertex_id = u + query_type_offset_;
        if (similar_degree[u] >= p_mu)
        {
            core_num++;
            if (print_num++ > 100)
            {
                continue;
            }
            cout << vertex_id << "-" << similar_degree[u] << "-" << qn_adj_List[u].size();
            cout << " visited: " << visited_qtv_[u] << endl;
            int sim_nei_num = 0;
            for (int i = 0; i < qn_adj_List[u].size(); i++)
            {
                int v = qn_adj_List[u][i] - query_type_offset_;
                bool sim_res = check_struc_sim(u, v);

                if (sim_res)
                {
                    sim_nei_num++;
                    // cout << qn_adj_List[u][i] << "-" << sim_res << " ";
                }
                else
                {
                    // cout << "n" << v << "-" << 0 << " ";
                }
            }
            // cout << endl;
            // cout << "--real sim num of QN: " << sim_nei_num << endl;
            if (sim_nei_num < p_mu)
            {
                cout << "-------------------------this vertex is not an core vertex---------------" << endl;
            }
        }
    }
    cout << "all query type vertices number is " << num_query_type_ << endl;
    cout << "the num of core vertices: " << core_num << endl;
    return;
}

void HinGraph::output_result(string output)
{
    string output_file_path = output + "/query_type_" + to_string(p_query_type) + "-mu" + to_string(p_mu);
    output_file_path += "-" + data_file_name + "-" + query_file_name + "-mode" + to_string(mode_query) + ".txt";
    ofstream output_file = open_file_ofstream(output_file_path);
    cout << "output file name and path: " << output_file_path << endl;
    output_file << "c/n vertex_id -- cluster_id\n";

    for (int i = 0; i < num_query_type_; i++)
    {
        int vertex_id = i + query_type_offset_;
        if (is_in_community[i])
        {
            if (cand_core_[i])
            {
                output_file << "c " << vertex_id << endl;
            }
            else
            {
                output_file << "n " << vertex_id << endl;
            }
        }
    }

    output_file.close();
    // fclose(fout);
}

int HinGraph::get_vertex_type(int vertex_id)
{
    if (vertex_id < 0 || vertex_id >= n)
    {
        printf("Can not get the vertex type: %d\n", vertex_id);
        exit(1);
    }
    int type = -1;
    for (int i = 0; i < n_types - 1; i++)
    {
        if (vertex_id >= vertex_start_map_[i] && vertex_id < vertex_start_map_[i + 1])
        {
            type = i;
            break;
        }
    }
    if (type == -1)
    {
        type = n_types - 1;
    }
    return type;
}