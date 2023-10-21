#include "HinGraph.h"

int initial_vector_size = 32;

bool judge_demoinate(const vector<float> &vec1, const vector<float> &vec2);

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
    data_index_dir_ = data_dir + "/index-cs";
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
    hin_schema_edge_cnt.resize(n_types, ::vector<int>(n_types, 0));
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
        mode_query = 10;
        if (mode == "-qidx1")
        {
            mode_query = 11;
        }
        if (mode == "-qidxON")
        {
            mode_query = 111;
        }
        index_query_();
        return;
    }
    else
    {
        if (mode == "-q")
        {
            baseline_query_();
            return;
        }

        if (mode == "-q1")
        {
            mode_query = 1;
        }
        else if (mode == "-q2")
        {
            mode_query = 2;
            baseline_query_();
            return;
        }
        else if (mode == "-q3")
        {
            mode_query = 3;
        }
        improved_query_();
    }

    // count_core_vertices();
}

void HinGraph::construct_index(string query_file, string option)
{
    k_max = 10;
    load_query_file(query_file);
    check_dir_path(data_index_dir_);
    initialize_query_();
    cout << "start construct index for vertex type: " << p_query_type << endl;
    Timer t1;
    t1.Start();

    search_d_neighbor();
    check_empty_set();
    save_d_n();
    t1.StopAndPrint("search k-strata time");
    // compute_all_similarity();
    // save_all_similarity();
    load_all_similarity();
    t1.StopAndPrint("compute similarity time");
    compute_k_threshold();
    return;
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
            query_node_list.resize(query_node_num);
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
            if (vertex_type == p_query_type)
            {
                distance_[vertex_type] = -1;
                continue;
            }
            type_epsilon[vertex_type] = unit_epsilon_value;
            cout << vertex_type << " -> " << type_epsilon[vertex_type] << " -- ";
            if (distance >= 1)
            {
                int type_start = vertex_start_map_[vertex_type];
                int type_end = (vertex_type + 1) == n_types ? n : vertex_start_map_[vertex_type + 1];
                // for (int j = type_start; j < type_end; j++)
                // {
                //     t_hop_visited[j] = vector<int>();
                //     t_hop_visited.reserve(initial_vector_size);
                // }
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
                // for (int j = type_start; j < type_end; j++)
                // { // create a visited set for each two neighbor
                //     t_hop_visited[j] = vector<int>();
                //     t_hop_visited.reserve(initial_vector_size);
                // }
                // cout << "create a visited unordered_set for each two neighbor :";
                // cout << type_start << "-" << type_end << endl;
            }
        }
    }
    int min_distance = 1000;
    for (auto &pair : type_epsilon)
    {
        if (pair.first == p_query_type || pair.second == 0.0)
        {
            auto it = type_epsilon.find(pair.first);
            if (it != type_epsilon.end())
                type_epsilon.erase(it);
            continue;
        }
        if (distance_[pair.first] <= min_distance)
        {
            cand_gen_type = pair.first;
            min_distance = distance_[pair.first];
        }
    }
    // cout << cand_gen_type << " ";

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
    res_community.resize(num_query_type_);
    similar_degree.resize(num_query_type_);
    effective_degree.resize(num_query_type_);
    cand_core_.resize(num_query_type_);
    node_k_thres.resize(num_query_type_);
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
        dn_adj_List[i].type_order_ = vector<int>();
        dn_adj_List[i].type_order_.reserve(n_types);
        qn_adj_List[i] = vector<int>();
        qn_adj_List[i].reserve(initial_vector_size);
        cand_sn_list[i] = vector<int>();
        cand_sn_list[i].reserve(initial_vector_size);
        non_sn_list[i] = vector<int>();
        non_sn_list[i].reserve(initial_vector_size);
        visited_qtv_[i] = false;
        ks_visit[i] = false;
        is_in_community[i] = false;
        res_community[i] = false;
        similar_degree[i] = 0;
        effective_degree[i] = 0;
        cand_core_[i] = false;
    }
    mlq.initial(num_query_type_);
    // cout << "finish reinitialize variable" << endl;
    return;
}

void HinGraph::print_result(bool print_all, long use_time)
{
    if (mode_query == 10)
    {
        if (is_in_community[query_i])
        {
            int num_com = 0;
            for (int i = 0; i < num_query_type_; i++)
            {
                if (is_in_community[i])
                {
                    num_com++;
                    if (print_all)
                        cout << i;
                }
            }
            cout << query_i + query_type_offset_ << " community contains Num(V) : " << num_com;
            cout << ". Use time: " << use_time << endl;
        }
        return;
    }

    if (similar_degree[query_i] >= p_mu)
    {
        if (mode_query == 1)
        {
            int num_community = 0, core_num = 0;
            for (int i = 0; i < num_query_type_; i++)
            {
                if (is_in_community[i])
                {
                    num_community++;
                    if (print_all)
                        cout << i;
                    if (cand_core_[i])
                    {
                        core_num++;
                        if (print_all)
                            cout << "c";
                    }
                    if (print_all)
                        cout << " ";
                }
            }

            cout << query_i + query_type_offset_ << " community contains Num(V) : " << num_community;
            cout << " Core number: " << core_num;
            cout << ". Use time: " << use_time << endl;
            return;
        }
        int num_community = 0;
        for (int i = 0; i < num_query_type_; i++)
        {
            if (res_community[i])
            {
                num_community++;
                if (print_all)
                    cout << i<<" ";
            }
        }

        cout << query_i + query_type_offset_ << " community contains Num(V) : " << num_community;
        cout << ". Use time: " << use_time << endl;
    }
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
        queue<int> cs_node_q, delete_q;
        is_in_community[query_i] = true;
        if (query_i + query_type_offset_ == 4860)
            cout << "error" << endl;

        Timer t1;
        t1.Start();
        search_k_strata(query_i);
        search_cand_sn(query_i);
        check_sn(query_i, cs_node_q, delete_q);
        if (similar_degree[query_i] < p_mu)
        {
            cout << query_i << " Cannot search a community" << endl;
        }
        else
        {
            while (!cs_node_q.empty())
            { // get Homo Graph
                int vertex_i = cs_node_q.front();
                cs_node_q.pop();
                search_cand_sn(vertex_i);
                check_sn(vertex_i, cs_node_q, delete_q);
            }
            core_decomposition(delete_q);
        }

        long cost_time = t1.StopTime();
        all_time += cost_time;
        time_cost[i] = cost_time;
        print_result(query_i + query_type_offset_ == 4860, cost_time);
        if (query_i + query_type_offset_ == 4860)
            break;
    }
    string str1 = "finish query " + to_string(query_node_num) + " times, use time:";
    Timer::PrintTime(str1, all_time);
}

void HinGraph::index_query_()
{
    load_all_similarity();
    load_k_thres_vec(p_mu);
    index_order_epsilon.resize(index_type_order.size());
    for (int i = 0; i < index_type_order.size(); i++)
    {
        int type_i = index_type_order[i];
        if (type_epsilon.find(type_i) == type_epsilon.end())
        {
            index_order_epsilon[i] = 0.0;
            continue;
        }
        if (type_epsilon[type_i] == 0.0)
        {
            index_order_epsilon[i] = 0.0;
            continue;
        }
        index_order_epsilon[i] = type_epsilon[type_i];
    }

    cout << "start index query " << endl;
    long all_time = 0;
    vector<long> time_cost(query_node_num, 0);
    for (int i = 0; i < query_node_num; i++)
    {
        query_i = query_node_list[i];
        for (int i = 0; i < query_node_num; i++)
            is_in_community[query_i] = false;

        queue<int> cs_node_q;
        Timer t1;
        t1.Start();

        if (index_judge_core(query_i, p_mu) == false)
        {
            cout << query_i << " Cannot search a community" << endl;
        }
        else
        {
            cs_node_q.push(query_i);
            while (!cs_node_q.empty())
            { // get all k-core
                int vertex_i = cs_node_q.front();
                cs_node_q.pop();
                for (const auto &i_nei : h_sim[vertex_i])
                {
                    if (is_in_community[i_nei.neighbor_i])
                        continue; // this neighbor has been add to community
                    if (judge_demoinate(i_nei.sim_vec, index_order_epsilon) == false)
                        continue; // this neighbor cannot add to current community
                    if (index_judge_core(i_nei.neighbor_i, p_mu) == false)
                        continue; // this neighbor is not a core
                    cs_node_q.push(i_nei.neighbor_i);
                    is_in_community[i_nei.neighbor_i] = true;
                }
            }
        }

        long cost_time = t1.StopTime();
        all_time += cost_time;
        time_cost[i] = cost_time;
        print_result(false, cost_time);
    }
    string str1 = "finish query " + to_string(query_node_num) + " times, use time:";
    Timer::PrintTime(str1, all_time);
}

bool HinGraph::index_judge_core(int i, int k)
{
    const auto &ik_thres = node_k_thres[i][k].thres_vecs;
    for (const auto &t_v : ik_thres)
    {
        if (judge_demoinate(t_v, index_order_epsilon))
        {
            is_in_community[i] = true;
            return true;
        }
    }
    is_in_community[i] = false;
    return false;
}

void HinGraph::estimate_best_order(vector<int> &order_type)
{
    queue<int> q;
    vector<vector<int>> bin_types(n_types, vector<int>());
    for (auto pair : type_epsilon)
    {
        if (pair.second == 0)
            continue;
        if (pair.first == p_query_type)
        {
            bin_types[2].push_back(pair.first);
            continue;
        }
        bin_types[distance_[pair.first]].push_back(pair.first);
    }
    vector<int> types_node_cnt(n_types, 0);
    for (int i = 0; i < n_types; i++)
    { // get the node count of each type
        int end_num = (i + 1) == n_types ? n : vertex_start_map_[i + 1];
        types_node_cnt[i] = end_num - vertex_start_map_[i];
    }
    vector<double> types_edge_ratio(n_types, 0);
    for (int i = 0; i < n_types; i++)
    {
        int cur_dis = i;
        if (bin_types[i].size() == 0)
            continue;
        if (bin_types[i].size() == 1)
        {
            int cur_type = bin_types[i][0];
            order_type.push_back(cur_type);
            if (cur_dis == 1)
            {
                types_edge_ratio[cur_type] = hin_schema_edge_cnt[cur_type][p_query_type];
            }
            else
            {
                for (int j = 0; j < n_types; j++)
                {
                    if (hin_schema_adjacencyMatrix[cur_type][j] == 0)
                        continue;
                    if (distance_[j] == cur_dis - 1)
                    {
                        double edge_num = hin_schema_edge_cnt[cur_type][j] * types_edge_ratio[j] / types_node_cnt[j];
                        types_edge_ratio[cur_type] += edge_num;
                    }
                }
            }
            continue;
        }

        // compute the order of same distance
        vector<pair<int, double>> same_dis_type_score;
        for (auto cur_type : bin_types[i])
        {
            if (cur_dis == 1)
            {
                types_edge_ratio[cur_type] = hin_schema_edge_cnt[cur_type][p_query_type];
            }
            else
            {
                for (int j = 0; j < n_types; j++)
                {
                    if (hin_schema_adjacencyMatrix[cur_type][j] == 0)
                        continue;
                    if (distance_[j] == cur_dis - 1)
                    {
                        double edge_num = hin_schema_edge_cnt[cur_type][j] * types_edge_ratio[j] / types_node_cnt[j];
                        types_edge_ratio[cur_type] += edge_num;
                    }
                }
            }
            same_dis_type_score.push_back(make_pair(cur_type, types_edge_ratio[cur_type] / types_node_cnt[cur_type]));
        }
        sort(same_dis_type_score.begin(), same_dis_type_score.end(), [](const auto &a, const auto &b)
             { return a.second < b.second; });
        for (auto pair : same_dis_type_score)
        {
            order_type.push_back(pair.first);
        }
    }
    cand_gen_type = order_type[0];
    for (auto type : order_type)
        cout << type << " ";
    cout << endl;
}

void HinGraph::improved_query_()
{
    vector<int> order_type;
    estimate_best_order(order_type);
    cout << "start improve online query " << endl;
    long all_time = 0;
    vector<long> time_cost(query_node_num, 0);
    for (int i = 0; i < query_node_num; i++)
    {
        query_i = query_node_list[i];
        reinitialize_query_();
        // queue<int> cs_node_q;
        is_in_community[query_i] = true;
        Timer t1;
        t1.Start();
        search_k_strata(query_i);
        search_cand_sn(query_i);
        cs_check_cluster_core_order(query_i, order_type);
        if (similar_degree[query_i] < p_mu)
        {
            cout << query_i << " Cannot search a community" << endl;
        }
        else
        {
            while (!mlq.empty())
            {
                int vertex_i = mlq.get_next();
                if (vertex_i == -1)
                    break;
                while (visited_qtv_[vertex_i]) // this vertex has been visited
                {
                    vertex_i = mlq.get_next();
                    if (vertex_i == -1)
                        break;
                }
                if (vertex_i == -1)
                    break;
                search_cand_sn(vertex_i);
                cs_check_cluster_core_order(vertex_i, order_type);
            }
        }

        long cost_time = t1.StopTime();
        all_time += cost_time;
        time_cost[i] = cost_time;
        print_result(false, cost_time);
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
    int can_gen = cand_gen_type;
    if (dn_i.num_d_neighbor[can_gen] == 0)
    {
        for (auto pair : type_epsilon)
        {
            if (pair.second == 0)
                continue;
            if (dn_i.num_d_neighbor[pair.first] != 0)
            {
                can_gen = pair.first;
                break;
            }
        }
    }
    const vector<int> &gen_type_dn = dn_i.d_neighbor_[can_gen];
    int cand_gen_dis = distance_[can_gen];
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
            // if (cur_step >= cand_gen_dis)
            if (cur_step >= p_d)
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
                    // if (distance_[nei_type] > cur_dis) // wrong direction
                    //     continue;
                    if (v_vertex.find(nei_id) != v_vertex.end()) // this neighbor visited before
                        continue;

                    v_vertex.insert(nei_id);
                    gen_q.push(make_pair(nei_id, cur_step + 1));
                }
            }
        }
    }
    sort(cand_sn_list[i].begin(), cand_sn_list[i].end());
    for (int cand_i : cand_sn_list[i])
    {
        int v = cand_i - query_type_offset_;
        if (ks_visit[v] == false)
        {
            search_k_strata(v);
        }
    }
}

void HinGraph::cs_check_cluster_core(int u, queue<int> &cs_queue)
{
    // 1. check core
    int vertex_u_id = u + query_type_offset_;
    int j = 0;
    unordered_set<int> v_vertex(qn_adj_List[u].begin(), qn_adj_List[u].end());
    for (auto not_i : non_sn_list[u])
        v_vertex.insert(not_i);

    effective_degree[u] = cand_sn_list[u].size();
    if (similar_degree[u] < p_mu)
    {
        similar_degree[u] = qn_adj_List[u].size();
        vector<int> tmp_in_com;
        tmp_in_com.reserve(cand_sn_list[u].size());
        for (; j < cand_sn_list[u].size(); j++)
        {
            int v_id = cand_sn_list[u][j];
            int v = v_id - query_type_offset_;
            if ((v_vertex.find(v_id) != v_vertex.end())) // this sn has been computed before
                continue;
            if ((mode_query == 2) && is_in_community[v] == true)
            {
                tmp_in_com.push_back(v_id);
                continue;
            }
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
        for (auto in_v_id : tmp_in_com)
        {
            int v = in_v_id - query_type_offset_;
            bool sim_res = check_struc_sim(u, v);
            if (sim_res)
            {
                similar_degree[u]++;
                qn_adj_List[u].push_back(in_v_id);
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

void HinGraph::check_sn(int u, queue<int> &cs_queue, queue<int> &delete_q)
{
    // 1. check core
    int vertex_u_id = u + query_type_offset_;
    unordered_set<int> v_vertex(qn_adj_List[u].begin(), qn_adj_List[u].end());
    for (auto not_i : non_sn_list[u])
        v_vertex.insert(not_i);

    similar_degree[u] = qn_adj_List[u].size();
    vector<int> tmp_in_com;
    tmp_in_com.reserve(cand_sn_list[u].size());
    for (int j = 0; j < cand_sn_list[u].size(); j++)
    {
        int v_id = cand_sn_list[u][j];
        int v = v_id - query_type_offset_;
        if ((v_vertex.find(v) != v_vertex.end())) // this sn has been computed before
            continue;
        bool sim_res = check_struc_sim(u, v);
        if (sim_res)
        {
            similar_degree[u]++;
            qn_adj_List[u].push_back(v);
        }
        if (visited_qtv_[v] == false && v != u)
        {
            if (sim_res)
            {
                similar_degree[v]++;
                qn_adj_List[v].push_back(u);
            }
            else
                non_sn_list[v].push_back(u);
        }
    }
    visited_qtv_[u] = true;

    // 2. cluster core
    if (similar_degree[u] < p_mu)
    { // this vertex neighbor cannot expanded
        cand_core_[u] = false;
        delete_q.push(u);
        return;
    }
    cand_core_[u] = true;
    for (auto sn_i : qn_adj_List[u])
    {
        if (is_in_community[sn_i] == false)
        {
            is_in_community[sn_i] = true; // add into community
            cs_queue.push(sn_i);          // add to expand queue
        }
    }
}

void HinGraph::core_decomposition(queue<int> &delete_q)
{
    bool community = true;
    while (!delete_q.empty())
    {
        int cur_u = delete_q.front();
        delete_q.pop();
        if (cur_u == query_i)
        {
            community = false;
            break;
        }
        if (is_in_community[cur_u] == false) // this vertex has been deleted before
            continue;
        for (const auto nei_u : qn_adj_List[cur_u])
        {
            if (similar_degree[nei_u] >= p_mu)
            {
                similar_degree[nei_u]--;
                if (similar_degree[nei_u] < p_mu)
                    delete_q.push(nei_u);
            }
        }
        is_in_community[cur_u] = false;
        similar_degree[cur_u] = -1;
    }
    if (community == false)
    {
        cout << query_i << " Cannot search a community" << endl;
        similar_degree[query_i] = -1;
        return;
    }
    queue<int> expand_res;
    expand_res.push(query_i);
    res_community[query_i] = true;
    while (!expand_res.empty())
    {
        int cur_u = expand_res.front();
        expand_res.pop();
        for (const auto nei_u : qn_adj_List[cur_u])
        {
            if (is_in_community[nei_u])
                if (res_community[nei_u] == false)
                { // connect k-core neighbor and not expand
                    res_community[nei_u] = true;
                    expand_res.push(nei_u);
                }
        }
    }
}

void HinGraph::cs_check_cluster_core_order(int u, vector<int> order_type)
{
    // 1. check core
    int vertex_u_id = u + query_type_offset_;
    int j = 0;
    unordered_set<int> v_vertex(qn_adj_List[u].begin(), qn_adj_List[u].end());
    for (auto not_i : non_sn_list[u])
        v_vertex.insert(not_i);

    // vector<pair<int, int>> type_num_;
    // for (const auto &pair : type_epsilon)
    // {
    //     if (pair.second == 0.0)
    //         continue;
    //     type_num_.push_back(make_pair(pair.first, dn_adj_List[u].num_d_neighbor[pair.first]));
    // }
    // sort(type_num_.begin(), type_num_.end(), [](const auto &a, const auto &b)
    //      { return a.second < b.second; });
    // for (auto pair : type_num_)
    //     dn_adj_List[u].type_order_.push_back(pair.first);

    if (similar_degree[u] < p_mu)
    {
        similar_degree[u] = qn_adj_List[u].size();
        effective_degree[u] = cand_sn_list[u].size();
        vector<int> tmp_in_com;
        tmp_in_com.reserve(cand_sn_list[u].size());
        for (; j < cand_sn_list[u].size(); j++)
        {
            int v_id = cand_sn_list[u][j];
            int v = v_id - query_type_offset_;
            if (v_vertex.find(v_id) != v_vertex.end()) // this sn has been computed before
                continue;
            if ((mode_query == 2) && is_in_community[v] == true)
            {
                tmp_in_com.push_back(v_id);
                continue;
            }

            // bool sim_res = check_struc_sim_with_order(u, v, dn_adj_List[u].type_order_);
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
                    // mlq.push(v, similar_degree[v]);
                }
                else
                    non_sn_list[v].push_back(vertex_u_id);
            }
            if (effective_degree[u] < p_mu || similar_degree[u] >= p_mu)
                break;
        }
        for (auto in_v_id : tmp_in_com)
        {
            int v = in_v_id - query_type_offset_;
            // bool sim_res = check_struc_sim_with_order(u, v, dn_adj_List[u].type_order_);
            bool sim_res = check_struc_sim(u, v);
            if (sim_res)
            {
                similar_degree[u]++;
                qn_adj_List[u].push_back(in_v_id);
            }
            else
                effective_degree[u]--;

            if (visited_qtv_[v] == false && v != u)
            {
                if (sim_res)
                {
                    similar_degree[v]++;
                    qn_adj_List[v].push_back(vertex_u_id);
                    // mlq.push(v, similar_degree[v]);
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
            // cs_queue.push(sn_i);          // add to expand queue
            mlq.push(sn_i, similar_degree[sn_i]);
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

        // bool sim_res = check_struc_sim_with_order(u, v, dn_adj_List[u].type_order_);
        bool sim_res = check_struc_sim(u, v);
        if (sim_res)
        {
            similar_degree[u]++;
            qn_adj_List[u].push_back(v_id);
        }

        if (visited_qtv_[v] == false && v != u)
        {
            if (sim_res)
            {
                similar_degree[v]++;
                qn_adj_List[v].push_back(vertex_u_id);
                if (is_in_community[v] == false)
                {
                    is_in_community[v] = true;
                    // cs_queue.push(v); // add to expand queue
                    mlq.push(v, similar_degree[v]);
                }
            }
            else
                non_sn_list[v].push_back(vertex_u_id);
        }
    }
    return;
    /*
    // 1. check core with order type
    int vertex_u_id = u + query_type_offset_;
    unordered_set<int> v_vertex(qn_adj_List[u].begin(), qn_adj_List[u].end());
    for (auto not_i : non_sn_list[u])
        v_vertex.insert(not_i);

    similar_degree[u] = qn_adj_List[u].size();
    vector<int> tmp;
    tmp.reserve(cand_sn_list[u].size());
    for (auto cand_i : cand_sn_list[u])
    {
        if (v_vertex.find(cand_i - query_type_offset_) != v_vertex.end()) // this sn has been computed before
            continue;
        tmp.push_back(cand_i);
    }
    cand_sn_list[u] = move(tmp);
    effective_degree[u] = cand_sn_list[u].size() + qn_adj_List[u].size();
    if (effective_degree[u] < p_mu)
    {
        cand_core_[u] = false;
        return;
    }
    int order_type_num = order_type.size();
    for (int i = 0; i < order_type_num; i++)
    {
        int type_i = order_type[i];
        vector<int> tmp;
        tmp.reserve(cand_sn_list[u].size());
        for (int j = 0; j < cand_sn_list[u].size(); j++)
        {
            int v_id = cand_sn_list[u][j];
            int v = v_id - query_type_offset_;
            bool sim_res = check_one_type_sim(u, v, type_i);

            if (sim_res)
                tmp.push_back(v_id);
            else
                effective_degree[u]--;

            if (i == (order_type_num - 1))
            {
                if (sim_res)
                {
                    similar_degree[u]++;
                    qn_adj_List[u].push_back(v_id);
                    if (visited_qtv_[v] == false && v != u)
                    {
                        similar_degree[v]++;
                        qn_adj_List[v].push_back(vertex_u_id);
                    }
                    if (is_in_community[v] == false)
                    {
                        is_in_community[v] = true;
                        cs_queue.push(v); // add to expand queue
                    }
                }
                else
                {
                    if (visited_qtv_[v] == false && v != u)
                        non_sn_list[v].push_back(vertex_u_id);
                }
            }

            if (effective_degree[u] < p_mu)
            {
                cand_core_[u] = false;
                return;
            }
        }
        cand_sn_list[u] = move(tmp);
    }

    if (similar_degree[u] < p_mu)
    {
        cand_core_[u] = false;
        return;
    }
    cand_core_[u] = true;
    */
}

void HinGraph::search_d_neighbor()
{
    cout << "search_all_d_neighborhood" << endl;
    cout << "all query type vertices number is " << num_query_type_ << endl;

    vector<int> type_degree_max(5, 0);
    vector<int> type_degree_min(5, 10);

    for (int i = 0; i < num_query_type_; i++)
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
                if (nei_dis >= 1)
                {
                    t_hop_visited[nei_id].push_back(query_vertex_id); // add to visited
                }
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
            if (dn_adj_List[i].d_neighbor_[type].size() == 0)
            {
                empty_dn_set[type].push_back(query_vertex_id);
            }
            else
            {
                sort(dn_adj_List[i].d_neighbor_[type].begin(), dn_adj_List[i].d_neighbor_[type].end());
                // edge_num[type] += dn_adj_List[i].num_d_neighbor[type];
                if (dn_adj_List[i].num_d_neighbor[type] > type_degree_max[type])
                    type_degree_max[type] = dn_adj_List[i].num_d_neighbor[type];
                if (dn_adj_List[i].num_d_neighbor[type] < type_degree_min[type])
                    type_degree_min[type] = dn_adj_List[i].num_d_neighbor[type];
            }
        }
        if (i % 100000 == 0)
            cout << i << " -- ";
    }
    cout << endl;
    for (auto &t_d : type_degree_max)
        cout << t_d << " ";
    cout << endl;
    for (auto &t_d : type_degree_min)
        cout << t_d << " ";
    cout << endl;
}

void HinGraph::save_d_n()
{
    cout << "start save d-neighbor" << endl;
    string type_index_path = data_index_dir_ + "/" + to_string(p_query_type);
    check_dir_path(type_index_path);
    string d_neighbor_path = type_index_path + "/d_neighbor1.txt";
    string d_neighbor2_path = type_index_path + "/d_neighbor2.txt";

    // 1. d_neighbor
    ofstream d_neighbor_file = open_file_ofstream(d_neighbor_path);
    for (int i = 0; i < num_query_type_; i++)
    {
        d_neighbor_file << i << " " << type_epsilon.size() << endl;
        for (const auto &pair : type_epsilon)
        {
            int type = pair.first;
            d_neighbor_file << type << " " << dn_adj_List[i].d_neighbor_[type].size() << " ";
            for (const auto &d_neighbor : dn_adj_List[i].d_neighbor_[type])
            {
                d_neighbor_file << d_neighbor << " ";
            }
            d_neighbor_file << endl;
        }
    }
    // 1.2 empty_dn_set
    d_neighbor_file << empty_dn_set.size() << endl;
    for (const auto &empty : empty_dn_set)
    {
        d_neighbor_file << empty.first << " " << empty.second.size() << endl;
        for (const auto &empty_x : empty.second)
        {
            d_neighbor_file << empty_x << " ";
        }
    }
    d_neighbor_file.close();

    // 2. d_neighbor-side adj-List
    ofstream d_neighbor2_file = open_file_ofstream(d_neighbor2_path);
    d_neighbor2_file << t_hop_visited.size() << endl;
    for (const auto &t_hop : t_hop_visited)
    {
        const int t_hop_vertex = t_hop.first;
        const vector<int> &query_type_n_list = t_hop.second;
        d_neighbor2_file << t_hop_vertex << " " << query_type_n_list.size() << endl;
        for (const auto &query_type_n : query_type_n_list)
        {
            d_neighbor2_file << query_type_n << " ";
        }
        d_neighbor2_file << endl;
    }
    d_neighbor2_file.close();

    return;
}

void HinGraph::load_d_n()
{
    cout << "start load d-neighbor" << endl;
    string type_index_path = data_index_dir_ + "/" + to_string(p_query_type);
    string d_neighbor_path = type_index_path + "/d_neighbor1.txt";
    string d_neighbor2_path = type_index_path + "/d_neighbor2.txt";

    // 1.1 d_neighbor
    ifstream d_neighbor_file = open_file_fstream(d_neighbor_path);
    for (int i = 0; i < num_query_type_; i++)
    {
        int load_i = -1, type_num = -1;
        d_neighbor_file >> load_i >> type_num;
        if (load_i != i)
        {
            cout << "load d neighbor file error" << endl;
            exit(-1);
        }
        for (int type_j = 0; type_j < type_num; type_j++)
        {
            int type, type_dn_len;
            d_neighbor_file >> type >> type_dn_len;
            vector<int> tmp(type_dn_len);
            for (int j = 0; j < type_dn_len; j++)
                d_neighbor_file >> tmp[j];
            if (type_epsilon.find(type) != type_epsilon.end())
            { // this type is considered in this query
                dn_adj_List[i].d_neighbor_[type] = move(tmp);
                dn_adj_List[i].num_d_neighbor[type] = type_dn_len;
            }
        }
    }
    // 1.2 empty_dn_set load
    int empty_size = 0;
    d_neighbor_file >> empty_size;
    for (int i = 0; i < empty_size; i++)
    {
        int empty_type_j, empty_type_j_size;
        d_neighbor_file >> empty_type_j >> empty_type_j_size;
        empty_dn_set[empty_type_j].reserve(empty_type_j_size);
        for (int j = 0; j < empty_type_j_size; j++)
        {
            int x;
            d_neighbor_file >> x;
            empty_dn_set[empty_type_j].push_back(x);
        }
    }
    d_neighbor_file.close();

    // 2. d_neighbor-side adj-List
    ifstream d_neighbor2_file = open_file_fstream(d_neighbor2_path);
    int d2_len;
    d_neighbor2_file >> d2_len;
    for (int i = 0; i < d2_len; i++)
    {
        int t_hop_vertex_id, t_hop_len;
        d_neighbor2_file >> t_hop_vertex_id >> t_hop_len;
        vector<int> tmp(t_hop_len);
        for (int j = 0; j < t_hop_len; j++)
        {
            d_neighbor2_file >> tmp[j];
        }
        if (t_hop_visited.find(t_hop_vertex_id) != t_hop_visited.end())
        { // this vertex is a t_hop_visited in this query
            t_hop_visited[t_hop_vertex_id] = move(tmp);
        }
    }
    d_neighbor2_file.close();

    return;
}

void HinGraph::compute_all_similarity()
{
    cout << "compute all similarity, with the type order: ";
    for (const auto pair : type_epsilon)
    {
        if (pair.second == 0)
            continue;
        index_type_order.push_back(pair.first);
        cout << pair.first << " ";
    }
    cout << endl;
    h_sim.resize(num_query_type_);

    for (int i = 0; i < num_query_type_; i++)
    {
        vector<Nei_similarity> qn_sim;
        for (int j = 0; j < index_type_order.size(); j++)
        {
            int type_j = index_type_order[j];
            vector<Query_nei_similarity> type_j_qn_similarity;
            compute_one_type_qn_similarity(i, type_j, type_j_qn_similarity);
            if (j == 0)
            {
                qn_sim.reserve(type_j_qn_similarity.size());
                for (const auto nei_sim : type_j_qn_similarity)
                    qn_sim.emplace_back(Nei_similarity{nei_sim.neighbor_i, nei_sim.similarity});
            }
            else
                intersection_neisim(qn_sim, type_j_qn_similarity);
        }
        compute_domin_rank(qn_sim);
        h_sim[i] = move(qn_sim);
        if (i % (num_query_type_ / 10) == 0)
        {
            cout << i << endl;
        }
    }
}

void HinGraph::save_all_similarity()
{
    cout << "start save all similarity graph" << endl;
    string type_index_path = data_index_dir_ + "/" + to_string(p_query_type);
    check_dir_path(type_index_path);
    string all_sim_path = type_index_path + "/sim_graph.txt";

    // 1. index_type_order
    ofstream all_sim_file = open_file_ofstream(all_sim_path);
    all_sim_file << num_query_type_ << " " << index_type_order.size() << endl;
    for (auto type : index_type_order)
    {
        all_sim_file << type << " ";
    }
    all_sim_file << endl;
    // 2. similarity homo graph
    for (int i = 0; i < num_query_type_; i++)
    {
        all_sim_file << i << " " << h_sim[i].size() << endl;
        for (const auto nei_sim_i : h_sim[i])
        {
            all_sim_file << nei_sim_i.neighbor_i << " " << nei_sim_i.domin_rank << " ";
            for (auto type_sim : nei_sim_i.sim_vec)
                all_sim_file << type_sim << " ";
            all_sim_file << endl;
        }
    }
    all_sim_file.close();
    cout << "finish save all similarity graph" << endl;
    return;
}

void HinGraph::load_all_similarity()
{
    string type_index_path = data_index_dir_ + "/" + to_string(p_query_type);
    string all_sim_path = type_index_path + "/sim_graph.txt";
    ifstream all_sim_file = open_file_fstream(all_sim_path);

    // 1. index_type_order
    int type_order_size;
    all_sim_file >> num_query_type_ >> type_order_size;
    index_type_order.resize(type_order_size);
    for (int i = 0; i < type_order_size; i++)
    {
        all_sim_file >> index_type_order[i];
    }
    // 2. similarity homo graph
    h_sim.resize(num_query_type_);
    for (int i = 0; i < num_query_type_; i++)
    {
        int input_i, i_nei_size;
        all_sim_file >> input_i >> i_nei_size;
        if (input_i != i)
        {
            cout << "load all similarity homo graph error" << endl;
            exit(-1);
        }

        h_sim[i].resize(i_nei_size);
        for (int j = 0; j < i_nei_size; j++)
        {
            all_sim_file >> h_sim[i][j].neighbor_i >> h_sim[i][j].domin_rank;
            h_sim[i][j].sim_vec.resize(type_order_size);
            for (int type_i = 0; type_i < type_order_size; type_i++)
                all_sim_file >> h_sim[i][j].sim_vec[type_i];
        }
    }
    all_sim_file.close();
    cout << "finish load all similarity graph" << endl;
    last_vertex.resize(num_query_type_);
    for (int i = 0; i < num_query_type_; i++)
        last_vertex[i] = true;
}

bool judge_demoinate(const vector<float> &vec1, const vector<float> &vec2)
{ // return vec1 dominate vec2 or not
    for (int i = 0; i < vec1.size(); i++)
    {
        if (vec1[i] < vec2[i])
            return false;
    }
    return true;
}

void HinGraph::compute_domin_rank(vector<Nei_similarity> &qn_sim)
{
    for (int i = 0; i < qn_sim.size(); i++)
    {
        int domin_cnt = 0;
        for (int j = 0; j < qn_sim.size(); j++)
        {
            if (i == j)
                continue;
            if (judge_demoinate(qn_sim[j].sim_vec, qn_sim[i].sim_vec))
                domin_cnt++;
        }
        qn_sim[i].domin_rank = domin_cnt + 1;
    }
}

bool generate_next_fix_thres(vector<float> &fix_type_thresold, vector<int> &fix_)
{ // remain_type is 0, and first type is 1 in fix_type_threshold
    fix_[0] += 1;
    for (int i = 0; i < fix_.size(); i++)
    {
        if ((fix_[i] == 100) && i != fix_.size() - 1)
        {
            fix_[i] = 0;
            fix_[i + 1] += 1;
        }
        if (i == fix_.size() - 1 && fix_[i] == 100)
            return true;
    }
    for (int i = 0; i < fix_.size(); i++)
    {
        fix_type_thresold[i] = fix_[i] / 100.0;
    }
    return false;
}

bool customSort(const std::vector<float> &a, const std::vector<float> &b)
{
    for (int dim = 0; dim < a.size(); ++dim)
    {
        if (a[dim] < b[dim])
        {
            return true;
        }
        else if (a[dim] > b[dim])
        {
            return false;
        }
    }
    return false; // a and b are equal
}

void sortVectorsByDimensions(std::vector<std::vector<float>> &thres_vecs)
{
    std::sort(thres_vecs.begin(), thres_vecs.end(), customSort);
}

void HinGraph::compute_k_threshold()
{
    cout << "start compute k threshold" << endl;
    node_k_thres.resize(num_query_type_);
    // p_mu = 1, 2 it contain itself and dominated value.
    for (int i = 0; i < num_query_type_; i++)
    {
        node_k_thres[i].reserve(h_sim[i].size());
        for (int k = 1; k <= 2; k++)
        {
            if (h_sim[i].size() < k)
                continue;
            if (k == 1 || k == 2)
            {
                k_threshold tmp;
                for (int j = 0; j < h_sim[i].size(); j++)
                {
                    if (h_sim[i][j].domin_rank == k)
                        tmp.thres_vecs.push_back(h_sim[i][j].sim_vec);
                }
                node_k_thres[i][k] = move(tmp);
            }
        }
    }

    // compute each_threshold
    for (int k = 3; k <= k_max; k++)
    {
        int re_type_offset = index_type_order.size() - 1;
        vector<float> fix_type_thresold(index_type_order.size() - 1, 0.0);
        vector<int> fix_(index_type_order.size() - 1, 0);
        load_all_similarity();
        cout << "current fix type and threshold:" << endl;
        for (int i = 0; i < fix_type_thresold.size(); i++)
            cout << i << "-" << fix_type_thresold[i] << "  ";
        cout << "|";
        compute_connect_k_core(k, fix_type_thresold, re_type_offset);
        while (generate_next_fix_thres(fix_type_thresold, fix_) == false)
        {
            bool reload_graph = fix_type_thresold[0] == 0.0 ? true : false;
            if (reload_graph == true)
            {
                load_all_similarity();
                cout << "current fix type and threshold is " << endl;
            }
            for (int i = 0; i < fix_type_thresold.size(); i++)
                cout << i << "-" << fix_type_thresold[i] << " ";
            cout << "|";
            if (fix_[0] % 10 == 0)
                cout << endl;
            compute_connect_k_core(k, fix_type_thresold, re_type_offset);
        }
        // for (int i = 0; i < num_query_type_; i++)
        // {
        //     sortVectorsByDimensions(node_k_thres[i][k].thres_vecs);
        // }
        save_k_thres_vec(k);
    }
}

bool judge_fix_type_edge(const vector<float> &sim_vec, const vector<float> &fix_type)
{
    for (int i = 0; i < fix_type.size(); i++)
    {
        if (sim_vec[i] < fix_type[i])
            return false;
    }

    return true;
}

void HinGraph::compute_connect_k_core(int k, const vector<float> &fix_type, int re_type)
{
    queue<int> delete_q;
    community_vertex.resize(num_query_type_, false);

    // 1. compute core number
    vector<int> coreNum(num_query_type_, 0); // core Number
    for (int i = 0; i < num_query_type_; i++)
    {
        if (last_vertex[i] == false) // this vertex has not been a core in current graph
            continue;

        vector<Nei_similarity> tmp_sim;
        tmp_sim.reserve(h_sim[i].size());
        for (const auto &i_nei_sim : h_sim[i])
        {
            if (judge_fix_type_edge(i_nei_sim.sim_vec, fix_type) == false)
            { // current edge sim not meet fix type threshold
                continue;
            }
            tmp_sim.push_back(i_nei_sim);
        }
        h_sim[i] = move(tmp_sim);
        if (h_sim[i].size() < k)
        { // current vertex cannot be a core at current threshold
            community_vertex[i] = false;
            delete_q.push(i);
            last_vertex[i] = false;
        }
        coreNum[i] = h_sim[i].size();
        community_vertex[i] = true; // current vertex maybe a initial k-core
        last_vertex[i] = true;
    }

    // 2. core decomposition
    while (!delete_q.empty())
    {
        int cur_i = delete_q.front();
        delete_q.pop();
        community_vertex[cur_i] = false;
        last_vertex[cur_i] = false;
        for (auto nei_cur_i : h_sim[cur_i])
        {
            int nei_i = nei_cur_i.neighbor_i;
            if (coreNum[nei_cur_i.neighbor_i] >= k)
            {
                coreNum[nei_cur_i.neighbor_i]--;
                if (coreNum[nei_cur_i.neighbor_i] < k)
                    delete_q.push(nei_cur_i.neighbor_i);
            }
        }
    }

    // 3. construct k-core graph for remain type change
    k_homo_graph.resize(num_query_type_, vector<k_homo_adj_node>());

    for (int i = 0; i < num_query_type_; i++)
    {
        if (community_vertex[i] == false)
            continue; // not core vertex
        k_homo_graph[i].clear();
        k_homo_graph[i].reserve(coreNum[i]);
        vector<Nei_similarity> tmp_sim;
        tmp_sim.reserve(coreNum[i]);
        for (const auto &i_nei_sim : h_sim[i])
        {
            if (community_vertex[i_nei_sim.neighbor_i] == true)
            {
                k_homo_graph[i].push_back(k_homo_adj_node(i_nei_sim.neighbor_i, i_nei_sim.sim_vec[re_type]));
                tmp_sim.push_back(i_nei_sim);
            }
        }
        h_sim[i] = move(tmp_sim);
        if (h_sim[i].size() != coreNum[i])
        {
            cout << i << " error " << coreNum[i] << " " << h_sim[i].size() << endl;
            exit(-1);
        }
    }

    // 4. compute all the remain type threshold
    for (int i = 0; i < 100; i++)
    {
        float re_type_threshold = i / 100.0;
        bool iter_next = search_and_add_threshold(k, fix_type, re_type_threshold);
        if (iter_next == false)
            break;
    }
    return;
}

void add_new_th(const vector<float> thres_sim_vec, vector<vector<float>> &old_thres_vecs)
{
    vector<vector<float>> tmp;
    tmp.reserve(old_thres_vecs.size());
    for (const auto old_vec : old_thres_vecs)
    {
        if (judge_demoinate(thres_sim_vec, old_vec))
            continue;
        tmp.push_back(old_vec); // old thres vec is not dominate by new one
    }
    tmp.push_back(thres_sim_vec);
    old_thres_vecs = move(tmp);
}

bool HinGraph::search_and_add_threshold(int k, const vector<float> &fix_type, float re_type_threshold)
{
    // 0. make similarity threshold vector
    vector<float> thres_sim_vec(fix_type);
    thres_sim_vec.push_back(re_type_threshold);

    // 1. compute core number
    queue<int> delete_q;
    vector<int> coreNum(num_query_type_, 0); // core Number
    for (int i = 0; i < num_query_type_; i++)
    {
        if (community_vertex[i] == false) // this vertex has not in current graph
            continue;
        vector<k_homo_adj_node> tmp_adj;
        tmp_adj.reserve(k_homo_graph[i].size());
        for (auto nei_j : k_homo_graph[i])
        {
            if (nei_j.re_type_sim > re_type_threshold) // only retain edge great than threshold. and cureent threshold is new thres
                tmp_adj.push_back(nei_j);
        }
        k_homo_graph[i] = move(tmp_adj);
        if (k_homo_graph[i].size() < k)
        { // current vertex cannot be a core at current threshold
            community_vertex[i] = false;
            delete_q.push(i);
        }
        coreNum[i] = k_homo_graph[i].size();
        community_vertex[i] = true; // current vertex maybe a initial k-core
    }

    // 2. core decomposition
    while (!delete_q.empty())
    {
        int cur_i = delete_q.front();
        delete_q.pop();
        community_vertex[cur_i] = false;
        add_new_th(thres_sim_vec, node_k_thres[cur_i][k].thres_vecs);
        for (auto nei_cur_i : k_homo_graph[cur_i])
        {
            if (coreNum[nei_cur_i.neighbor_i] >= k)
            {
                coreNum[nei_cur_i.neighbor_i]--;
                if (coreNum[nei_cur_i.neighbor_i] < k)
                    delete_q.push(nei_cur_i.neighbor_i);
            }
        }
    }

    // 3. construct k-core graph for next iteration
    bool next_iteration = false;
    for (int i = 0; i < num_query_type_; i++)
    {
        if (community_vertex[i] == false)
            continue; // not core vertex
        vector<k_homo_adj_node> tmp_adj;
        tmp_adj.reserve(coreNum[i]);
        for (auto &nei : k_homo_graph[i])
        {
            if (community_vertex[nei.neighbor_i] == true)
                tmp_adj.push_back(nei);
        }
        if (tmp_adj.size() != coreNum[i])
        {
            cout << " error " << endl;
            exit(-1);
        }
        k_homo_graph[i] = move(tmp_adj);
        next_iteration = true;
    }
    return next_iteration;
}

void HinGraph::intersection_neisim(vector<Nei_similarity> &nei_sim, const vector<Query_nei_similarity> &vec2)
{
    vector<Nei_similarity> intersection;
    intersection.reserve(nei_sim.size());
    int i = 0; // 指向 vec1 的指针
    int j = 0; // 指向 vec2 的指针

    while (i < nei_sim.size() && j < vec2.size())
    {
        if (nei_sim[i].neighbor_i == vec2[j].neighbor_i)
        {
            nei_sim[i].sim_vec.push_back(vec2[j].similarity);
            intersection.push_back(nei_sim[i]);
            i++;
            j++;
        }
        else if (nei_sim[i].neighbor_i < vec2[j].neighbor_i)
        {
            i++;
        }
        else
        {
            j++;
        }
    }
    nei_sim = move(intersection);
    return;
}

int HinGraph::compute_one_type_qn_similarity(int i, int type_i, vector<Query_nei_similarity> &type_i_qn_similarity)
{ // for construct index
    vector<int> &dn_i_type_i = dn_adj_List[i].d_neighbor_[type_i];
    if (dn_i_type_i.size() == 0)
    {
        if (mode_query == 0)
            return 0;
        type_i_qn_similarity.reserve(empty_dn_set[type_i].size());
        for (const auto &qn_i : empty_dn_set[type_i])
        {
            type_i_qn_similarity.emplace_back(Query_nei_similarity{qn_i - query_type_offset_, 1.0f});
        }
        return empty_dn_set[type_i].size();
    }

    vector<int> mergedVector;
    bool query_qn = (distance_[type_i] == 0);
    concat_one_type_qn(dn_i_type_i, mergedVector, query_qn);

    // computer the similarity
    type_i_qn_similarity.reserve(mergedVector.size());
    int prev_id = mergedVector[0], count_id_cnt = 1, all_type_i_qn = 0;
    int len_i = dn_adj_List[i].d_neighbor_[type_i].size();
    for (int j = 1; j < mergedVector.size(); j++)
    {
        if (prev_id == mergedVector[j])
        {
            count_id_cnt++;
        }
        else
        {
            all_type_i_qn++;
            int len_j = dn_adj_List[prev_id - query_type_offset_].d_neighbor_[type_i].size();
            double sim = static_cast<double>(count_id_cnt) / (len_i + len_j - count_id_cnt);
            type_i_qn_similarity.emplace_back(Query_nei_similarity{prev_id - query_type_offset_, static_cast<float>(sim)});

            prev_id = mergedVector[j];
            count_id_cnt = 1;
        }
    }
    int len_j = dn_adj_List[prev_id - query_type_offset_].d_neighbor_[type_i].size();
    double sim = static_cast<double>(count_id_cnt) / (len_i + len_j - count_id_cnt);
    all_type_i_qn++;
    type_i_qn_similarity.emplace_back(Query_nei_similarity{prev_id - query_type_offset_, static_cast<float>(sim)});
    return all_type_i_qn;
}

void HinGraph::concat_one_type_qn(const vector<int> &dn_i_type_i, vector<int> &mVector, bool query_type)
{
    size_t totalSize = 0;
    vector<int> mergedVector;

    for (const auto &d_n_i : dn_i_type_i)
    {
        const vector<int> &induced_qn = query_type ? dn_adj_List[d_n_i - query_type_offset_].d_neighbor_[query_type] : t_hop_visited[d_n_i];
        totalSize += induced_qn.size();
    }
    mergedVector.reserve(totalSize);
    for (const auto &d_n_i : dn_i_type_i)
    {
        const vector<int> &induced_qn = query_type ? dn_adj_List[d_n_i - query_type_offset_].d_neighbor_[query_type] : t_hop_visited[d_n_i];
        // const vector<int> &induced_qn = t_hop_visited[d_n_i];
        for (const auto &induced_qn_1 : induced_qn)
        {
            mergedVector.push_back(induced_qn_1);
        }
        // mergedVector.insert(mergedVector.end(), induced_qn.begin(), induced_qn.end());
    }
    sort(mergedVector.begin(), mergedVector.end());
    mVector = move(mergedVector);
}

void HinGraph::save_k_thres_vec(int k)
{
    cout << "start save threshold vector " << k << endl;
    string type_index_path = data_index_dir_ + "/" + to_string(p_query_type);
    check_dir_path(type_index_path);
    string k_thres_path = type_index_path + "/k_thres_" + to_string(k) + ".txt";

    ofstream k_thres_file = open_file_ofstream(k_thres_path);
    k_thres_file << num_query_type_ << " " << k << " " << index_type_order.size() << endl;
    for (int i = 0; i < num_query_type_; i++)
    {
        k_thres_file << i << " " << node_k_thres[i][k].thres_vecs.size() << endl;
        for (const auto &thres_vec : node_k_thres[i][k].thres_vecs)
        {
            for (const auto &vec_i : thres_vec)
            {
                k_thres_file << vec_i << " ";
            }
            k_thres_file << endl;
        }
    }
    cout << "finish save threshold vector " << k << endl;
}

void HinGraph::load_k_thres_vec(int k)
{
    string type_index_path = data_index_dir_ + "/" + to_string(p_query_type);
    string k_thres_path = type_index_path + "/k_thres_" + to_string(k) + ".txt";

    ifstream k_thres_file = open_file_fstream(k_thres_path);
    if (!k_thres_file)
    {
        cerr << "Error: Failed to open k_thres file for reading." << endl;
        return;
    }

    int num_query_types, k_value, num_index_types;
    k_thres_file >> num_query_types >> k_value >> num_index_types;
    node_k_thres.resize(num_query_type_);

    if (num_query_types != num_query_type_ || k_value != k || num_index_types != index_type_order.size())
    {
        cerr << "Error: Incompatible data found in k_thres file." << endl;
        return;
    }
    for (int i = 0; i < num_query_type_; i++)
    {
        int query_type, num_i_tv;
        k_thres_file >> query_type >> num_i_tv;

        if (query_type != i)
        {
            cerr << "Error: Incompatible data found for query type " << i << " in k_thres file." << endl;
            return;
        }
        node_k_thres[i][k].thres_vecs.resize(num_i_tv);
        for (int j = 0; j < num_i_tv; j++)
        {
            node_k_thres[i][k].thres_vecs[j].resize(num_index_types);
            for (int vec_i = 0; vec_i < num_index_types; vec_i++)
            {
                k_thres_file >> node_k_thres[i][k].thres_vecs[j][vec_i];
            }
        }
    }
}

bool HinGraph::check_struc_sim(int a, int b)
{
    if (a == b)
        return true;

    for (const auto &pair : type_epsilon)
    {
        int type = pair.first;
        double type_ep = pair.second;
        if (type_ep == 0.0)
            continue;

        if (judgeJacSim(dn_adj_List[a].d_neighbor_[type], dn_adj_List[b].d_neighbor_[type], type_ep) == false)
            return false;
    }
    return true;
}

bool HinGraph::check_struc_sim_with_order(int a, int b, const vector<int> order_type)
{
    if (a == b)
        return true;

    for (const auto &type : order_type)
    {
        double type_ep = type_epsilon[type];
        if (type_ep == 0.0)
            continue;

        if (judgeJacSim(dn_adj_List[a].d_neighbor_[type], dn_adj_List[b].d_neighbor_[type], type_ep) == false)
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

bool HinGraph::check_one_type_sim(int a, int b, int type)
{
    if (a == b)
        return true;
    double type_ep = type_epsilon[type];
    if (type_ep == 0)
        return true;
    if (judgeJacSim(dn_adj_List[a].d_neighbor_[type], dn_adj_List[b].d_neighbor_[type], type_ep))
        return true;
    else
        return false;
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

void HinGraph::check_empty_set()
{
    for (const auto &pair : type_epsilon)
    {
        int type = pair.first;
        double type_epsilon = pair.second;
        if (type_epsilon == 0.0 || distance_[type] < 1)
        {
            continue;
        }
        set<int> uniqueValues(empty_dn_set[type].begin(), empty_dn_set[type].end());

        if (empty_dn_set[type].size() != uniqueValues.size())
        {
            cout << "type " << type << " empty query vertices contains duplicate numbers. " << endl;
        }
        cout << "type " << type << " empty query vertices number is " << empty_dn_set[type].size() << endl;
    }
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