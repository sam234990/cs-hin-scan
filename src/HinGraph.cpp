#include "HinGraph.h"

using namespace std;

int initial_vector_size = 32;

bool judge_demoinate(const vector<float> &vec1, const vector<float> &vec2);

double local_clustering_coefficient(const vector<vector<int>> &adj_list, const vector<int> &res, int node);
double global_clustering_coefficient(const vector<vector<int>> &adj_list, const vector<int> &res);

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
    vec1 = std::move(intersection);
    return;
}

vector<int> randomSampling(int n, int sampleSize)
{
    vector<int> samples;
    int interval;
    if (n <= sampleSize)
    {
        interval = 1;
        cout << "the vertex number of current type is less than " << sampleSize;
        cout << ", use all the vertices" << endl;
    }
    else
        interval = n / sampleSize; // internal
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
    input_dir_ = data_dir;
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
        iss_graph_type >> edges_[i].v_id >> edges_[i].v_type >> edges_[i].edge_type;
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

void HinGraph::cs_hin_scan(string query_file, string mode, int scale)
{
    option_ = mode;
    load_query_file(query_file);
    initialize_query_();
    reinitialize_query_();

    if (mode == "-qidxscal" || mode == "-qscal")
    {
        scale_ = scale;
        max_scale_num = num_query_type_ * scale / 100;
        num_query_type_ = max_scale_num;
        max_scale_id = max_scale_num + query_type_offset_;
        // query_node_list = move(randomSampling(num_query_type_, query_node_num));
    }

    if (mode == "-qidx" || mode == "-qidxsel" || mode == "-qidxSCAN" || mode == "-qidxscal" || mode == "-qidxeff")
    {

        mode_query = 12;
        if (mode == "-qidxsel")
        {
            mode_query = 11;
            index_cd();
            return;
        }
        if (mode == "-qidxSCAN")
        {
            mode_query = 12;
        }
        if (mode == "-qidxscal")
        {
            mode_query = 16;
        }
        if (mode == "-qidxscal")
        {
            mode_query = 15;
        }
        index_query_();
        return;
    }
    else
    {
        if (mode == "-q")
        {
            mode_query = 2;
            baseline_query_();
            return;
        }

        if (mode == "-qpath_sim")
        {
            mode_query = 1;
            baseline_pathsim_query_();
            return;
        }
        else if (mode == "-qSCAN")
        {
            mode_query = 2;
            baseline_query_();
            return;
        }
        else if (mode == "-qcdSCAN")
        {
            mode_query = 2;
            online_cd();
            return;
        }
        else if (mode == "-qscal")
        {
            mode_query = 6;
            baseline_query_();
            return;
        }
        else if (mode == "-qeff" || mode == "-qeffcos" || mode == "-qeffpathsim")
        {
            mode_query = 2;
            baseline_query_();
            return;
        }
    }
}

void HinGraph::construct_index(string query_file, string option, int start_k, int scale)
{
    option_ = option;
    k_max = start_k;
    load_query_file(query_file);
    check_dir_path(data_index_dir_);
    initialize_query_();
    initial_construct_index();
    cout << "start construct index for vertex type: " << p_query_type << endl;
    Timer t1;
    t1.Start();

    if (option == "-fidxscal")
    {
        // cout << "scale test:" << scale << endl;
        scale_ = scale;
        max_scale_num = num_query_type_ * scale / 100;
        num_query_type_ = max_scale_num;
        max_scale_id = max_scale_num + query_type_offset_;
    }

    string type_index_path = data_index_dir_ + "/" + to_string(p_query_type);
    string all_sim_path = type_index_path + "/sim_graph.txt";
    std::ifstream sim_file(all_sim_path.c_str());
    if (option == "-fidxscal" || option == "-fidx1scal" || option == "-fidx_DBP" || sim_file.good() == false)
    {
        cout << "start counput the Similarity graph" << endl;
        search_d_neighbor();
        check_empty_set();
        if (option == "-fidx_DBP")
            return;
        save_d_n();
        t1.StopAndPrint("search k-strata time");
        compute_all_similarity();
        save_all_similarity();
    }
    else
    {
        cout << "Already compute the Similarity graph" << endl;
    }

    load_all_similarity();
    t1.StopAndPrint("compute similarity time");

    if (option == "-fidx" || option == "-fidxscal")
    {
        cout << "used basic index construction" << endl;
        compute_k_threshold(start_k);
    }
    else if (option == "-fidx1" || option == "-fidx1scal")
    {
        mode_query = 101;
        cout << "used advanced index construction" << endl;
        improve_k_thres(start_k);
    }
    else if (option == "-fidxmutree")
    {
        mode_query = 105;
        build_tree_index(start_k);
    }

    return;
}

void HinGraph::find_meta(int type)
{
    cout << "Start find meta-path" << endl;
    p_query_type = type;
    // initial query neighbor
    query_type_offset_ = vertex_start_map_[p_query_type];
    int end = ((p_query_type + 1) == n_types) ? n : vertex_start_map_[p_query_type + 1];
    num_query_type_ = end - query_type_offset_;
    cout << num_query_type_ << endl;

    vector<long> res(1000000, 0);
    for (int i = 0; i < num_query_type_; i++)
    {
        int query_vertex_id = i + query_type_offset_;
        int start = vertex_offset_[query_vertex_id];
        int end = ((query_vertex_id + 1) == n) ? m : vertex_offset_[query_vertex_id + 1];
        for (int j = start; j < end; j++)
        {
            int nei_type = edges_[j].v_type, nei_id = edges_[j].v_id;
            int nei_edge_start = vertex_offset_[nei_id];
            int nei_end = ((nei_id + 1) == n) ? m : vertex_offset_[nei_id + 1];
            for (int k = nei_edge_start; k < nei_end; k++)
            {
                int type_2 = edges_[k].v_type;
                if (type_2 == nei_type || type_2 == p_query_type)
                    continue;
                int add_meta = (nei_type + 1) * 1000 + type_2 + 1;
                res[add_meta] += 1;
            }
        }

        if (i % (num_query_type_ / 100) == 0)
            cout << i << endl;
    }
    vector<int> indices(1000000);
    for (int i = 0; i < 1000000; ++i)
        indices[i] = i;
    sort(indices.begin(), indices.end(), [&res](int a, int b)
         { return res[a] > res[b]; });
    int print_i = 0;
    for (int i : indices)
    {
        if (res[i] == 0)
            continue;
        int type_1 = i / 1000 - 1;
        int type_2 = i % 1000 - 1;
        std::cout << type_1 << "\t" << type_2 << "\t" << res[i] << std::endl;
        print_i++;
        if (print_i > 200)
            break;
    }
    cout << "start save meta-path" << endl;
    string save_dir = "xxx";
    check_dir_path(save_dir);
    string save_path = save_dir + "/" + to_string(p_query_type) + ".txt";
    ofstream save_file = open_file_ofstream(save_path);
    print_i = 0;
    for (int i : indices)
    {
        if (res[i] == 0)
            continue;
        int type_1 = i / 1000 - 1;
        int type_2 = i % 1000 - 1;
        save_file << type_1 << "\t" << type_2 << "\t" << res[i] << endl;
        print_i++;
        if (print_i > 200)
            break;
    }
    cout << "finish save meta-path" << endl;
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
    bool selected_ids = false;
    string line;
    ifstream query_file = open_file_fstream(query_file_path);
    int lineCount = 0;
    int key_number = 0;
    random_query = false;
    query_pathsim = false;
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
        { // 3. query_type & query times / number of query nodes
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
                { // use random ids
                    random_query = true;
                    cout << "use random query id" << endl;
                    break;
                }
                if (i == 0 && query_node_list[i] == -2)
                { // use selected ids
                    selected_ids = true;
                    cout << "use selected query id" << endl;
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
                else if (key == -2 && value > 0)
                {
                    unit_epsilon = true;
                    unit_epsilon_value = value;
                    query_pathsim = true;
                    break;
                }
                else if (key == -3 && value > 0)
                {
                    unit_epsilon = true;
                    unit_epsilon_value = value;
                    query_pathsim = true;
                    cout << "Use the unit epsilon vector and eval Pathsim" << endl;
                    break;
                }
                else if (key == -4 && value > 0)
                {
                    query_pathsim = true;
                    key_number = int(value);
                    cout << "Use the special epsilon vector and eval Pathsim" << endl;
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
    if (key_number > 0)
    {
        string key_line;
        for (int i = 0; i < key_number; i++)
        {
            getline(query_file, key_line);
            istringstream iss2(key_line);
            int key;
            double value;
            iss2 >> key >> value;
            type_epsilon[key] = value;
        }
    }

    if (query_pathsim == true)
    {
        int path_num;
        string path_num_line;
        getline(query_file, path_num_line);
        path_num = stoi(path_num_line);
        for (int i = 0; i < path_num; ++i)
        {
            string metapath;
            float epsilon;
            if (getline(query_file, metapath))
            {
                string line;
                getline(query_file, line);
                epsilon = stof(line);
                if (epsilon == 0.0)
                {
                    std::cerr << "Failed to read epsilon from file." << std::endl;
                    exit(1);
                }
            }
            else
            {
                std::cerr << "Failed to read a group from file." << std::endl;
                exit(1);
            }

            metapath_vecs.push_back(metapath);
            pathsim_epsilon.push_back(epsilon);
            cout << epsilon << " --";
        }
        cout << endl;
    }

    query_file.close();

    if (selected_ids)
    { // use selected ids
        string selected_ids_path;
        if (option_ == "-qidxscal" || option_ == "-qscal")
            selected_ids_path = data_index_dir_ + "/" + to_string(p_query_type) + "/selected_scale_ids.txt";
        else
            selected_ids_path = data_index_dir_ + "/" + to_string(p_query_type) + "/selected_ids.txt";
        cout << "use Selected node for query " << selected_ids_path << endl;
        ifstream selected_ids_file = open_file_fstream(selected_ids_path);
        for (int i = 0; i < query_node_num; i++)
        {
            selected_ids_file >> query_node_list[i];
        }
    }

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
        std::vector<int> keysToRemove;
        for (auto it = type_epsilon.begin(); it != type_epsilon.end(); ++it)
        {
            if (it->second == 0.0)
                keysToRemove.push_back(it->first);
        }

        for (int key : keysToRemove)
            type_epsilon.erase(key);

        for (const auto &pair : type_epsilon)
        {
            cout << pair.first << " -> " << pair.second << " -- ";
        }
        cout << endl;
    }

    query_type_offset_ = vertex_start_map_[p_query_type];
    int end = ((p_query_type + 1) == n_types) ? n : vertex_start_map_[p_query_type + 1];
    num_query_type_ = end - query_type_offset_;

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
        }
        cout << endl;
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

    if (random_query)
    {
        query_node_list = move(randomSampling(num_query_type_, query_node_num));
    }

    if (random_query == false)
    {
        for (int i = 0; i < query_node_num && i < query_node_list.size(); i++)
        {
            if (query_node_list[i] > query_type_offset_)
                query_node_list[i] -= query_type_offset_;
        }
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
    similar_degree.resize(num_query_type_, 0);
    effective_degree.resize(num_query_type_);
    cand_core_.resize(num_query_type_);
    cout << "finish initial" << endl;
    return;
}

void HinGraph::initial_construct_index()
{
    for (const auto &entry : type_epsilon)
    {
        // int vertex_type = entry.first, distance = entry.second;
        int vertex_type = entry.first;
        // if (distance >= 1)
        int type_start = vertex_start_map_[vertex_type];
        int type_end = (vertex_type + 1) == n_types ? n : vertex_start_map_[vertex_type + 1];
        for (int j = type_start; j < type_end; j++)
        {
            int nei_edge_start = vertex_offset_[j];
            int nei_end = ((j + 1) == n) ? m : vertex_offset_[j + 1];
            t_hop_visited[j] = vector<int>();
            t_hop_visited[j].reserve((nei_end - nei_edge_start + 3));
        }
    }
    node_k_thres.resize(num_query_type_);
    community_vertex.resize(num_query_type_, false);
    for (int i = 0; i < num_query_type_; i++)
    {
        node_k_thres[i].thres_vecs.reserve(initial_vector_size);
        node_k_thres[i].corner_points.reserve(initial_vector_size);
    }
    cout << "neighbor type order: ";
    index_type_order.clear();
    for (const auto pair : type_epsilon)
    {
        if (pair.second == 0)
            continue;
        index_type_order.push_back(pair.first);
        empty_dn_set[pair.first].clear();
        empty_dn_set[pair.first].reserve(num_query_type_);
        cout << pair.first << " ";
    }
    cout << endl;
    index_tree.resize(index_type_order.size());
    k_homo_graph.resize(num_query_type_);
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
    if (mode_query >= 10)
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
                        cout << i << " ";
                }
            }
            if (print_all)
                cout << endl;
            cout << query_i << " community contains Num(V) : " << num_com;
            cout << ". Use time: " << use_time << endl;
        }
        return;
    }

    vector<int> res;
    if (is_in_community[query_i])
    {
        int num_community = 0;
        if (mode_query == 2 || mode_query == 1)
        {
            int core_num = 0;
            for (int i = 0; i < num_query_type_; i++)
            {
                if (is_in_community[i])
                {
                    num_community++;
                    if (query_node_num == 1)
                        res.push_back(i);
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

            cout << query_i << " community contains Num(V) : " << num_community;
            // cout << " Core number: " << core_num;
            cout << ". Use time: " << use_time << endl;
        }
        else
        {
            for (int i = 0; i < num_query_type_; i++)
            {
                if (res_community[i])
                {
                    num_community++;
                    if (print_all)
                        cout << i << " ";
                }
            }
            if (print_all)
                cout << endl;
            cout << query_i << " community contains Num(V) : " << num_community;
            cout << ". Use time: " << use_time << endl;
        }
    }
    if (query_node_num == 1)
    {
        vector<string> author_names(num_query_type_);
        string k_thres_path = data_dir_ + "/author.txt";
        ifstream inputFile(k_thres_path);
        std::string line;
        int idx = 0;
        while (std::getline(inputFile, line))
        {
            line.erase(std::remove(line.begin(), line.end(), '\n'), line.end());
            author_names[idx] = line;
            idx++;
        }
        for (auto i : res)
        {
            // cout << i << " ";
            cout << author_names[i] << ": ";
            for (auto j : res)
            {
                bool sim_res = check_struc_sim(i, j);
                if (sim_res)
                    cout << j << " ";
            }
            cout << endl;
        }
    }
}

void HinGraph::baseline_query_()
{
    cout << "start baseline online query " << endl;
    long all_time = 0;
    has_community = 0;
    vector<int> vertex_num_all(1000, 0);
    vector<int> core_num_all(1000, 0);
    vector<int> diameter_all(1000, 0);
    vector<double> density_all(1000, 0);
    vector<double> cc_all(1000, 0);
    vector<double> sim_all(1000, 0);
    vector<double> cos_all(1000, 0);
    vector<double> pathsim_all(1000, 0);
    if (option_ == "-qeffpathsim")
    {
        path_utils.initial_metapaths(metapath_vecs);
        path_utils.initial_query_vertex(num_query_type_);
    }

    vector<long> time_cost(query_node_num, 0);
    for (int i = 0; i < query_node_num && i < query_node_list.size(); i++)
    {
        query_i = query_node_list[i];
        reinitialize_query_();
        queue<int> cs_node_q, delete_q;
        is_in_community[query_i] = true;

        Timer t1;
        t1.Start();
        if (mode_query == 2)
        {
            online_query_scan();
        }
        else
        {
            search_k_strata(query_i);
            search_cand_sn(query_i);
            check_sn(query_i, cs_node_q, delete_q);
            if (similar_degree[query_i] < p_mu)
            {
                cout << query_i << " Cannot search a community" << endl;
                is_in_community[query_i] = false;
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
        }

        long cost_time = t1.StopTime();
        all_time += cost_time;
        time_cost[i] = cost_time;
        print_result(false, cost_time);
        if (option_ == "-qeff" || option_ == "-qeffcos" || option_ == "-qeffpathsim")
            online_effective_result(i, vertex_num_all, core_num_all, diameter_all,
                                    density_all, cc_all, sim_all, cos_all, pathsim_all);
    }
    string str1 = "finish query " + to_string(query_node_num) + " times, use time:";
    Timer::PrintTime(str1, all_time);
    if (option_ == "-qeff" || option_ == "-qeffcos" || option_ == "-qeffpathsim")
    {
        long vertex_num_ = 0, d_ = 0, core_num_ = 0;
        double den_ = 0.0, cc_ = 0.0, sim_ = 0.0, pathsim_ = 0.0;
        int com_num = 0;
        for (int i = 0; i < query_node_num && i < query_node_list.size(); i++)
        {
            if (vertex_num_all[i] == 0)
                continue;
            com_num++;
            vertex_num_ += vertex_num_all[i];
            core_num_ += core_num_all[i];
            d_ += diameter_all[i];
            den_ += density_all[i];
            cc_ += cc_all[i];
            sim_ += sim_all[i];
            pathsim_ += pathsim_all[i];
        }
        cout << "average vertex num is " << double(vertex_num_) / com_num << endl;
        cout << "average core num is " << double(core_num_) / com_num << endl;
        cout << "average diameter is " << double(d_) / com_num << endl;
        cout << "average density is " << den_ / com_num << endl;
        cout << "average cc is " << cc_ / com_num << endl;
        cout << "average similarity is " << sim_ / com_num << endl;
        cout << "average PathSim is " << pathsim_ / com_num << endl;
    }
}

void HinGraph::baseline_pathsim_query_()
{
    cout << "start baseline online query (use PathSim)" << endl;
    path_utils.initial_metapaths(metapath_vecs);
    path_utils.initial_query_vertex(num_query_type_);
    long all_time = 0;
    has_community = 0;
    vector<int> vertex_num_all(1000, 0);
    vector<int> core_num_all(1000, 0);
    vector<int> diameter_all(1000, 0);
    vector<double> density_all(1000, 0);
    vector<double> cc_all(1000, 0);
    vector<double> sim_all(1000, 0);
    vector<double> cos_all(1000, 0);
    vector<double> pathsim_all(1000, 0);

    vector<long> time_cost(query_node_num, 0);
    for (int i = 0; i < query_node_num && i < query_node_list.size(); i++)
    {
        query_i = query_node_list[i];
        // reinitialize_query_();
        for (int i = 0; i < num_query_type_; i++)
        {
            qn_adj_List[i] = vector<int>();
            qn_adj_List[i].reserve(initial_vector_size);
            cand_sn_list[i] = vector<int>();
            cand_sn_list[i].reserve(initial_vector_size);
            non_sn_list[i] = vector<int>();
            non_sn_list[i].reserve(initial_vector_size);
            visited_qtv_[i] = false;
            // ks_visit[i] = false;
            is_in_community[i] = false;
            res_community[i] = false;
            similar_degree[i] = 0;
            effective_degree[i] = 0;
            cand_core_[i] = false;
        }
        mlq.initial(num_query_type_);
        // path_utils.initial_query_vertex(num_query_type_);
        is_in_community[query_i] = true;
        Timer t1;
        t1.Start();

        path_utils.search(*this, query_i);
        path_utils.generate_cand_nei(*this, query_i, cand_sn_list[query_i]);
        scan_check_cluster_core(query_i);
        if (similar_degree[query_i] < p_mu)
        {
            is_in_community[query_i] = false;
            cout << query_i << " Cannot search a SCAN-like community (PathSim)" << endl;
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
                path_utils.generate_cand_nei(*this, vertex_i, cand_sn_list[vertex_i]);
                scan_check_cluster_core(vertex_i);
            }
        }
        long cost_time = t1.StopTime();
        all_time += cost_time;
        time_cost[i] = cost_time;
        print_result(false, cost_time);
        online_effective_result(i, vertex_num_all, core_num_all, diameter_all,
                                density_all, cc_all, sim_all, cos_all, pathsim_all);
    }
    int vertex_num_ = 0, d_ = 0, core_num_ = 0;
    double den_ = 0.0, cc_ = 0.0;
    double jacsim_ = 0.0, cossim_ = 0.0, pathsim_ = 0.0;
    int com_num = 0;
    for (int i = 0; i < query_node_num && i < query_node_list.size(); i++)
    {
        if (vertex_num_all[i] == 0)
            continue;
        com_num++;
        vertex_num_ += vertex_num_all[i];
        core_num_ += core_num_all[i];
        d_ += diameter_all[i];
        den_ += density_all[i];
        cc_ += cc_all[i];
        jacsim_ += sim_all[i];
        cossim_ += cos_all[i];
        pathsim_ += pathsim_all[i];
    }
    cout << "average vertex num is " << double(vertex_num_) / com_num << endl;
    cout << "average core num is " << double(core_num_) / com_num << endl;
    cout << "average diameter is " << double(d_) / com_num << endl;
    cout << "average density is " << den_ / com_num << endl;
    cout << "average cc is " << cc_ / com_num << endl;
    cout << "average jac similarity is " << jacsim_ / com_num << endl;
    cout << "average cos similarity is " << cossim_ / com_num << endl;
    cout << "average pathsim similarity is " << pathsim_ / com_num << endl;
}

void HinGraph::online_query_scan()
{
    if (ks_visit[query_i] == false)
    {
        search_k_strata(query_i);
    }
    search_cand_sn(query_i);
    scan_check_cluster_core(query_i);
    if (similar_degree[query_i] < p_mu)
    {
        is_in_community[query_i] = false;
    }
    else
    {
        cd_res_ssc.push_back(query_i);
        while (!mlq.empty())
        {
            int vertex_i = mlq.get_next();
            if (vertex_i == -1)
                break;
            cd_res_ssc.push_back(vertex_i);
            while (visited_qtv_[vertex_i]) // this vertex has been visited
            {
                vertex_i = mlq.get_next();
                if (vertex_i == -1)
                    break;
            }
            if (vertex_i == -1)
                break;
            search_cand_sn(vertex_i);
            scan_check_cluster_core(vertex_i);
        }
    }
}

void HinGraph::online_cd()
{
    cout << "start online community detection" << endl;
    long all_time = 0;
    has_community = 0;
    reinitialize_query_();
    cd_res_ssc.clear();
    cd_res_ssc.reserve(num_query_type_);
    vector<vector<int>> all_res_com;
    int community_all_num = 0;

    vector<int> c_member_i;
    c_member_i.reserve(2 * num_query_type_);
    vector<int> community_number(num_query_type_, 0);
    community_member.resize(num_query_type_);
    for (int i = 0; i < num_query_type_; i++)
        community_member[i] = false;

    Timer t0;
    t0.Start();
    vector<long> time_cost(num_query_type_, 0);
    for (int i = 0; i < num_query_type_; i++)
    {
        if (visited_qtv_[i] == true) // this vertex has been check core or not
            continue;
        query_i = i;

        is_in_community[query_i] = true;

        Timer t1;
        t1.Start();
        online_query_scan();

        long cost_time = t1.StopTime();
        all_time += cost_time;
        time_cost[i] = cost_time;
        print_result(false, cost_time);
        is_in_community[query_i] = false;

        if (cd_res_ssc.size() > 0)
        {
            int ssc_id = query_i;
            for (auto com_i : cd_res_ssc)
            {
                is_in_community[com_i] = false;
                community_number[com_i] = ssc_id;
                community_member[com_i] = true;
            }
            community_all_num++;
            sort(cd_res_ssc.begin(), cd_res_ssc.end());
            c_member_i.insert(c_member_i.end(), cd_res_ssc.begin(), cd_res_ssc.end());
            all_res_com.push_back(cd_res_ssc);
            mlq.initial(num_query_type_);
            cd_res_ssc.clear();
            cd_res_ssc.reserve(num_query_type_);
        }
    }

    cout << "finish find all SSC, start find hub and outlier" << endl;

    Others online_others;
    online_others.compute_hub_outlier_online(*this, c_member_i, community_number);

    cout << all_time << "contain community numbers:" << community_all_num << endl;

    int outlier_number = 0, hub_number = 0;
    for (int i = 0; i < num_query_type_; i++)
    {
        if (online_others.outlier[i] == true)
            outlier_number++;
        if (online_others.hub[i] == true)
            hub_number++;
    }
    outlier_number -= hub_number;
    cout << "outlier number: " << outlier_number << "  hub number :" << hub_number << endl;

    string str1 = "finish online community detection, use time:";
    t0.StopAndPrint(str1);
    // Timer::PrintTime(str1, all_time);
}

void HinGraph::scan_check_cluster_core(int u)
{
    // 1. check core
    int vertex_u_id = u + query_type_offset_;
    int j = 0;
    vector<bool> v_vertex(n, false);
    for (auto qi : qn_adj_List[u])
        v_vertex[qi] = true;

    for (auto not_i : non_sn_list[u])
        v_vertex[not_i] = true;

    if (similar_degree[u] < p_mu)
    {
        similar_degree[u] = qn_adj_List[u].size();
        effective_degree[u] = cand_sn_list[u].size();
        // vector<int> tmp_in_com;
        // tmp_in_com.reserve(cand_sn_list[u].size());
        for (; j < cand_sn_list[u].size(); j++)
        {
            int v_id = cand_sn_list[u][j];
            int v = v_id - query_type_offset_;
            if (v_vertex[v_id] == true)
                continue;
            // if (v_vertex.find(v_id) != v_vertex.end()) // this sn has been computed before
            //     continue;
            // if ((mode_query == 2) && is_in_community[v] == true)
            // {
            //     tmp_in_com.push_back(v_id);
            //     continue;
            // }
            bool sim_res;
            if (option_ == "-qpath_sim")
                sim_res = path_utils.judge_pathsim(u, v, pathsim_epsilon);
            else
                sim_res = check_struc_sim(u, v);

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
            mlq.push(sn_i, similar_degree[sn_i]);
        }
    }
    for (; j < cand_sn_list[u].size(); j++)
    {
        int v_id = cand_sn_list[u][j];
        int v = v_id - query_type_offset_;
        if (is_in_community[v] == true) // already in community
            continue;

        if (v_vertex[v_id] == true)
            continue;
        // if (v_vertex.find(v_id) != v_vertex.end()) // this sn has been computed before
        //     continue;
        // if (ks_visit[v] == false) // this vertex has not searched k-strata
        // search_k_strata(v);

        bool sim_res;
        if (mode_query == 1)
            sim_res = path_utils.judge_pathsim(u, v, pathsim_epsilon);
        else
            sim_res = check_struc_sim(u, v);

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
                    mlq.push(v, similar_degree[v]);
                }
            }
            else
                non_sn_list[v].push_back(vertex_u_id);
        }
    }
    return;
}

void HinGraph::query_index_scan()
{
    queue<int> cs_node_q;

    if (index_judge_core(query_i, p_mu) == false)
    {
        int nei_num = 0;
        for (const auto &i_nei : h_sim[query_i])
        {
            if (judge_demoinate(i_nei.sim_vec, index_order_epsilon) == false)
                continue; // this neighbor cannot add to current community
            nei_num++;
        }
        if (nei_num < p_mu)
        {
            cout << query_i << " Cannot search a SCAN-like community" << endl;
            return;
        }
    }
    cs_node_q.push(query_i);
    is_in_community[query_i] = true;

    while (!cs_node_q.empty())
    {
        int vertex_i = cs_node_q.front();
        cs_node_q.pop();
        if (index_judge_core(vertex_i, p_mu))
        {
            cand_core_[vertex_i] = true;
            for (const auto &i_nei : h_sim[vertex_i])
            {
                if (is_in_community[i_nei.neighbor_i])
                    continue; // this neighbor has been add to community
                if (judge_demoinate(i_nei.sim_vec, index_order_epsilon) == false)
                    continue; // this neighbor cannot add to current community
                cs_node_q.push(i_nei.neighbor_i);
                is_in_community[i_nei.neighbor_i] = true;
            }
        }
        else
        {
            vector<int> nei_ids;
            nei_ids.reserve(p_mu);
            for (const auto &i_nei : h_sim[vertex_i])
            {
                if (judge_demoinate(i_nei.sim_vec, index_order_epsilon) == false)
                    continue; // this neighbor cannot add to current community
                nei_ids.push_back(i_nei.neighbor_i);
            }
            if (nei_ids.size() >= p_mu)
            {
                cand_core_[vertex_i] = true;
                for (const auto nei_i : nei_ids)
                {
                    if (is_in_community[nei_i] == true)
                        continue;
                    cs_node_q.push(nei_i);
                    is_in_community[nei_i] = true;
                }
            }
        }
    }
    has_community++;
}

int compute_diameter(const vector<vector<int>> &adj, const vector<int> &res)
{
    int diameter = 0;
    int roundThreshold = 10;
    int round = 0;
    vector<int> dia_res = res;
    random_shuffle(dia_res.begin(), dia_res.end());

    unordered_set<int> res_set(res.begin(), res.end());

    for (const auto &com_i : dia_res)
    {
        int cur_max_hop = 0;
        unordered_set<int> v_vertex;
        queue<pair<int, int>> bfs_path;
        bfs_path.push(make_pair(com_i, 0));
        v_vertex.insert(com_i);

        while (!bfs_path.empty())
        {
            int curVertex_id = bfs_path.front().first;
            int cur_step = bfs_path.front().second;
            bfs_path.pop();

            if (cur_step > cur_max_hop)
                cur_max_hop = cur_step;

            for (const auto &nei_id : adj[curVertex_id])
            {
                if (v_vertex.find(nei_id) != v_vertex.end())
                    continue; // This neighbor has been visited before
                if (res_set.find(nei_id) == res_set.end())
                    continue; // This nei not in res

                v_vertex.insert(nei_id);
                bfs_path.push(make_pair(nei_id, cur_step + 1));
            }
        }

        if (cur_max_hop > diameter)
            diameter = cur_max_hop;

        round++;
        if (round >= roundThreshold)
            break;
    }

    return diameter;
}

void HinGraph::effectiveness_result(int eff_res_i, vector<int> &vertex_num_all,
                                    vector<int> &diameter_all, vector<double> &density_all,
                                    vector<int> &core_num_all, vector<double> &cc_all,
                                    vector<double> &jac_all, vector<double> &cosine_all,
                                    vector<double> &pathsim_all, vector<int> &eff_id)
{
    if (is_in_community[query_i] == false)
        return;
    int num_com = 0, edge_num = 0, core_num = 0;
    vector<int> res;
    res.reserve(num_query_type_);
    for (int i = 0; i < num_query_type_; i++)
    {
        if (is_in_community[i])
        {
            num_com++;
            res.push_back(i);
        }
    }

    // density
    vector<int> core_ids;
    core_ids.reserve(num_com);
    long p_edge_num = 0;
    for (const auto &com_i : res)
    {
        if (is_in_community[com_i] == false)
            continue;
        if (cand_core_[com_i])
        {
            core_ids.push_back(com_i);
            core_num++;
        }
        // qn_adj_List[com_i].reserve(h_sim[com_i].size());
        for (const auto &i_nei : h_sim[com_i])
        {
            if (is_in_community[i_nei.neighbor_i] == false)
                continue; // this neighbor not in community
            if (judge_demoinate(i_nei.sim_vec, index_order_epsilon) == false)
                continue; // this edge should be ignore
            edge_num++;
            // qn_adj_List[com_i].push_back(i_nei.neighbor_i);
        }

        if (qn_adj_List[com_i].size() == 0)
        {
            qn_adj_List[com_i] = path_utils.p_induced_graph(*this, com_i);
        }
    }
    double density = double(edge_num) / num_com;
    double pdensity = double(p_edge_num) / num_com;
    density_all[eff_res_i] = density;
    vertex_num_all[eff_res_i] = num_com;
    core_num_all[eff_res_i] = core_num;
    // pdensity_all[eff_res_i] = pdensity;

    // diameter
    int diameter = compute_diameter(qn_adj_List, res);
    diameter_all[eff_res_i] = diameter;

    // cc
    double global_cc = global_clustering_coefficient(qn_adj_List, res);
    cc_all[eff_res_i] = global_cc;

    for (const auto &com_i : res)
    {
        qn_adj_List[com_i].clear();
    }

    double jac = 0, cos = 0, pm = 0;
    int sim_cnt = 0, sample_cnt = 0;
    vector<int> sampled_nodes;
    int sample_size = 100;
    if (res.size() > sample_size)
    {
        vector<int> indices(res.size());
        iota(indices.begin(), indices.end(), 0);
        shuffle(indices.begin(), indices.end(), mt19937{random_device{}()});
        sampled_nodes.resize(sample_size);
        copy(res.begin(), res.begin() + sample_size, sampled_nodes.begin());
    }
    else
    {
        sampled_nodes = res;
    }
    for (const auto &res_i : sampled_nodes)
    {
        path_utils.search(*this, res_i);
        for (const auto &res_j : res)
        {
            path_utils.search(*this, res_j);
            pm += path_utils.compute_avg_pathsim(*this, res_i, res_j);
            vector<double> sim_jc = avgjac_cosSim(res_i, res_j);
            jac += sim_jc[0];
            cos += sim_jc[1];
            sim_cnt++;
        }
    }

    jac_all[eff_res_i] = jac / sim_cnt;
    cosine_all[eff_res_i] = cos / sim_cnt;
    pathsim_all[eff_res_i] = pm / sim_cnt;

    int id1 = query_i;
    int id2 = -1;
    for (const auto &i_nei : h_sim[id1])
    {
        if (is_in_community[i_nei.neighbor_i] == false)
            continue; // this neighbor not in community
        if (judge_demoinate(i_nei.sim_vec, index_order_epsilon) == false)
            continue; // this edge should be ignore
        if (i_nei.neighbor_i <= id1)
            continue;
        id2 = i_nei.neighbor_i;
        break;
    }
    eff_id.push_back(id1);
    if (id2 == -1)
    {
        random_shuffle(core_ids.begin(), core_ids.end());
        int id = query_i;
        auto it = std::find_if(core_ids.begin(), core_ids.end(), [id](int val)
                               { return val != id; });
        int core_2 = *it;
        eff_id.push_back(core_2);
    }
    else
    {
        eff_id.push_back(id2);
        return;
    }
}

void HinGraph::online_effective_result(int eff_res_i, vector<int> &vertex_num_all, vector<int> &core_num_all,
                                       vector<int> &diameter_all, vector<double> &density_all,
                                       vector<double> &cc_all, vector<double> &sim_all,
                                       vector<double> &cos_all, vector<double> &pathsim_all)
{
    if (is_in_community[query_i] == false)
        return;
    // 1. vertex/core num
    int num_com = 0, edge_num = 0, core_num = 0;
    vector<int> res;
    res.reserve(num_query_type_);
    for (int i = 0; i < num_query_type_; i++)
    {
        if (is_in_community[i])
        {
            num_com++;
            res.push_back(i);
            if (cand_core_[i])
                core_num++;
        }
    }
    vertex_num_all[eff_res_i] = num_com;
    core_num_all[eff_res_i] = core_num;

    // 2. density & similarity -- structural similar edge
    double sim = 0.0;
    for (const int &com_i : res)
    {
        if (is_in_community[com_i] == false)
            continue;
        search_k_strata(com_i);
        qn_adj_List[com_i].clear();
        qn_adj_List[com_i] = path_utils.p_induced_graph(*this, com_i);
    }
    cout << "finish induce homo graph" << endl;
    double density = double(edge_num) / num_com;
    density_all[eff_res_i] = density;

    // 3. diameter
    int diameter = compute_diameter(qn_adj_List, res);
    diameter_all[eff_res_i] = diameter;

    // 4. clustering coefficient
    double global_cc = global_clustering_coefficient(qn_adj_List, res);
    cc_all[eff_res_i] = global_cc;

    for (const auto &com_i : res)
    {
        qn_adj_List[com_i].clear();
    }

    double jac = 0, cos = 0, pm = 0;
    int sim_cnt = 0, sample_cnt = 0;
    vector<int> sampled_nodes;
    int sample_size = 100;
    if (res.size() > sample_size)
    {
        vector<int> indices(res.size());
        iota(indices.begin(), indices.end(), 0);
        shuffle(indices.begin(), indices.end(), mt19937{random_device{}()});
        sampled_nodes.resize(sample_size);
        copy(res.begin(), res.begin() + sample_size, sampled_nodes.begin());
    }
    else
    {
        sampled_nodes = res;
    }
    for (const auto &res_i : sampled_nodes)
    {
        path_utils.search(*this, res_i);
        for (const auto &res_j : res)
        {
            path_utils.search(*this, res_j);
            pm += path_utils.compute_avg_pathsim(*this, res_i, res_j);
            vector<double> sim_jc = avgjac_cosSim(res_i, res_j);
            jac += sim_jc[0];
            cos += sim_jc[1];
            sim_cnt++;
        }
    }

    sim_all[eff_res_i] = jac / sim_cnt;
    cos_all[eff_res_i] = cos / sim_cnt;
    pathsim_all[eff_res_i] = pm / sim_cnt;
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
        int type_start = vertex_start_map_[type_i];
        int type_end = (type_i + 1) == n_types ? n : vertex_start_map_[type_i + 1];
        for (int j = type_start; j < type_end; j++)
        {
            int nei_edge_start = vertex_offset_[j];
            int nei_end = ((j + 1) == n) ? m : vertex_offset_[j + 1];
            t_hop_visited[j] = vector<int>();
            t_hop_visited[j].reserve((nei_end - nei_edge_start + 3));
        }
    }
    load_d_n();

    path_utils.initial_metapaths(metapath_vecs);
    path_utils.initial_query_vertex(num_query_type_);

    cout << "start index query " << endl;
    long all_time = 0;
    has_community = 0;
    vector<int> vertex_num_all(1000, 0);
    vector<int> core_num_all(1000, 0);
    vector<int> diameter_all(1000, 0);
    vector<double> density_all(1000, 0);
    vector<double> pdensity_all(1000, 0);
    vector<double> cc_all(1000, 0);
    vector<double> jac_all(1000, 0);
    vector<double> cos_all(1000, 0);
    vector<double> pathsim_all(1000, 0);
    vector<int> eff_id;
    eff_id.reserve(1000);

    vector<long> time_cost(query_node_num, 0);
    for (int i = 0; i < query_node_num && i < query_node_list.size(); i++)
    {
        query_i = query_node_list[i];
        for (int i = 0; i < num_query_type_; i++)
        {
            cand_core_[i] = false;
            is_in_community[i] = false;
        }

        queue<int> cs_node_q;
        Timer t1;
        t1.Start();

        if (mode_query == 12 || mode_query == 16 || mode_query == 15)
        {
            query_index_scan();
        }
        else
        {
            if (index_judge_core(query_i, p_mu) == false)
            {
                cout << query_i << " Cannot search a community" << endl;
                is_in_community[query_i] = false;
            }
            else
            {
                cs_node_q.push(query_i);
                is_in_community[query_i] = true;
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
                has_community += 1;
            }
        }

        long cost_time = t1.StopTime();
        if (option_ == "-qidxeff")
            effectiveness_result(i, vertex_num_all, diameter_all, density_all,
                                 core_num_all, cc_all, jac_all, cos_all, pathsim_all, eff_id);
        all_time += cost_time;
        time_cost[i] = cost_time;
        print_result(false, cost_time);
    }
    string str1 = "finish query " + to_string(query_node_num) + " times, use time:";
    Timer::PrintTime(str1, all_time);
    cout << "the number of search community:" << has_community << endl;

    if (option_ == "qidx")
        return;

    long vertex_num_ = 0, d_ = 0, core_ = 0;
    double den_ = 0.0;
    double pden_ = 0.0;
    double cc_ = 0.0, jac_ = 0.0, cos_ = 0.0, pathsim_ = 0.0;
    int com_num = 0;
    for (int i = 0; i < query_node_num && i < query_node_list.size(); i++)
    {
        if (vertex_num_all[i] == 0)
            continue;
        com_num++;
        vertex_num_ += vertex_num_all[i];
        core_ += core_num_all[i];
        d_ += diameter_all[i];
        den_ += density_all[i];
        pden_ += pdensity_all[i];
        cc_ += cc_all[i];
        jac_ += jac_all[i];
        cos_ += cos_all[i];
        pathsim_ += pathsim_all[i];
    }
    double average_num = double(vertex_num_) / com_num;
    double average_dia = double(d_) / com_num;
    double average_den = den_ / com_num;
    double average_pden = pden_ / com_num;
    double average_cc = cc_ / com_num;
    cout << "average vertex num is " << average_num << endl;
    cout << "average core num is " << double(core_) / com_num << endl;
    cout << "average diameter is " << average_dia << endl;
    cout << "average density is " << average_den << endl;
    cout << "average pdensity is " << average_pden << endl;
    cout << "average cc is " << average_cc << endl;
    cout << "average jaccard is " << jac_ / com_num << endl;
    cout << "average cosine is " << cos_ / com_num << endl;
    cout << "average pathsim is " << pathsim_ / com_num << endl;
    for (size_t i = 0; i < eff_id.size(); i++)
    {
        cout << eff_id[i] + query_type_offset_;
        if (i % 2 == 1)
            cout << endl;
        else
            cout << " ";
    }
}

float sum_float_vec(const vector<float> &vec)
{
    float sum = 0;
    for (const auto v : vec)
        sum += v;
    return sum;
}

double local_clustering_coefficient(const vector<vector<int>> &adj_list, const unordered_set<int> &res, int node)
{
    vector<int> common_neighbors;
    common_neighbors.reserve(adj_list[node].size());
    for (int neighbor : adj_list[node])
    {
        if (res.find(neighbor) != res.end())
            common_neighbors.push_back(neighbor);
    }
    int neighbors = common_neighbors.size();
    if (neighbors < 2)
    {
        return 0.0;
    }

    int triangle_count = 0;
    for (int i = 0; i < neighbors; ++i)
    {
        int u = common_neighbors[i];
        for (int j = i + 1; j < neighbors; ++j)
        {
            int v = common_neighbors[j];
            if (find(adj_list[u].begin(), adj_list[u].end(), v) != adj_list[u].end())
            {
                triangle_count++;
            }
        }
    }

    int possible_triangles = neighbors * (neighbors - 1) / 2;
    if (possible_triangles == 0)
    {
        return 0.0;
    }

    return static_cast<double>(triangle_count) / possible_triangles;
}

double global_clustering_coefficient(const vector<vector<int>> &adj_list, const vector<int> &res)
{
    vector<int> sampled_nodes;
    int sample_size = 100;
    if (res.size() > sample_size)
    {
        vector<int> indices(res.size());
        iota(indices.begin(), indices.end(), 0);
        shuffle(indices.begin(), indices.end(), mt19937{random_device{}()});
        sampled_nodes.resize(sample_size);
        copy(res.begin(), res.begin() + sample_size, sampled_nodes.begin());
    }
    else
    {
        sampled_nodes = res;
    }

    unordered_set<int> res_set(res.begin(), res.end());

    double total_cc = 0.0;
    for (int node : sampled_nodes)
    {
        total_cc += local_clustering_coefficient(adj_list, res_set, node);
    }

    return total_cc / res.size();
}

void HinGraph::index_cd()
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

    path_utils.initial_metapaths(metapath_vecs);
    path_utils.initial_query_vertex(num_query_type_);

    vector<int> community_number(num_query_type_, 0);
    cand_core_.resize(num_query_type_);
    pa.resize(num_query_type_);
    p_rank_.resize(num_query_type_);
    community_member.resize(num_query_type_);
    qn_adj_List.resize(num_query_type_);
    int bs_max_ed = -1, bs_max_ed_i;
    unordered_map<int, queue<int>> bin_sort_queue_;
    for (int i = 0; i < num_query_type_; i++)
    {
        cand_core_[i] = false;
        community_member[i] = false;
        is_in_community[i] = false;
        pa[i] = i;
        p_rank_[i] = 0;
        similar_degree[i] = 0;
        effective_degree[i] = h_sim[i].size();
        visited_qtv_[i] = false;                     // mark each vertex unvisited
        qn_adj_List[i].reserve(initial_vector_size); // for effective analysis
        if (effective_degree[i] >= p_mu)
            bin_sort_queue_[effective_degree[i]].push(i);
        if (effective_degree[i] > bs_max_ed)
            bs_max_ed = effective_degree[i];
    }

    int community_all_num = 0;
    vector<int> com_size;
    com_size.reserve(num_query_type_);
    vector<float> sim_all;
    sim_all.reserve(num_query_type_);
    vector<int> edge_all;
    edge_all.reserve(num_query_type_);
    vector<int> c_member_i;
    c_member_i.reserve(2 * num_query_type_);


    Timer t1;
    t1.Start();
    while (true)
    {
        query_i = -1;
        if (bs_max_ed < p_mu)
            break;
        if (!bin_sort_queue_[bs_max_ed].empty())
        {
            int tmp = bin_sort_queue_[bs_max_ed].front();
            bin_sort_queue_[bs_max_ed].pop();
            if (visited_qtv_[tmp] == true)
                continue;
            else
                query_i = tmp;
        }
        else
        {
            bs_max_ed--;
            continue;
        }
        if (query_i == -1)
            break;

        visited_qtv_[query_i] = true;
        queue<int> cs_node_q;
        if (index_judge_core(query_i, p_mu) == false)
        {
            int nei_num = 0;
            for (const auto &i_nei : h_sim[query_i])
            {
                if (judge_demoinate(i_nei.sim_vec, index_order_epsilon) == false)
                    continue; // this neighbor cannot add to current community
                nei_num++;
            }
            if (nei_num < p_mu)
                continue;
        }
        vector<int> res;
        res.reserve(num_query_type_);
        cs_node_q.push(query_i);
        is_in_community[query_i] = true;
        cand_core_[query_i] = true;

        while (!cs_node_q.empty())
        {
            int vertex_i = cs_node_q.front();
            cs_node_q.pop();
            visited_qtv_[vertex_i] = true;
            res.push_back(vertex_i);

            if (index_judge_core(vertex_i, p_mu))
            {
                cand_core_[vertex_i] = true;
                for (const auto &i_nei : h_sim[vertex_i])
                {
                    // if (find_root(vertex_i) == find_root(i_nei.neighbor_i))
                    if (is_in_community[i_nei.neighbor_i])
                        continue; // this neighbor has been add to community
                    if (judge_demoinate(i_nei.sim_vec, index_order_epsilon) == false)
                        continue; // this neighbor cannot add to current community
                    cs_node_q.push(i_nei.neighbor_i);
                    is_in_community[i_nei.neighbor_i] = true;
                    // my_union(vertex_i, i_nei.neighbor_i);
                }
            }
            else
            {
                vector<int> nei_ids;
                nei_ids.reserve(p_mu);
                for (const auto &i_nei : h_sim[vertex_i])
                {
                    if (judge_demoinate(i_nei.sim_vec, index_order_epsilon) == false)
                        continue; // this neighbor cannot add to current community
                    nei_ids.push_back(i_nei.neighbor_i);
                }
                if (nei_ids.size() >= p_mu)
                {
                    cand_core_[vertex_i] = true;
                    for (const auto nei_i : nei_ids)
                    {
                        // if (find_root(vertex_i) == find_root(nei_i))
                        if (is_in_community[nei_i] == true)
                            continue;
                        cs_node_q.push(nei_i);
                        is_in_community[nei_i] = true;
                        // my_union(vertex_i, nei_i);
                    }
                }
            }
        }
        int com_num = res.size();

        float community_sim = 0;
        int edge_num = 0;
        for (const auto &com_i : res)
        {
            if (cand_core_[com_i] == false)
                continue;
            for (const auto &i_nei : h_sim[com_i])
            {
                if (is_in_community[i_nei.neighbor_i] == false)
                    continue; // this neighbor not in community
                if (judge_demoinate(i_nei.sim_vec, index_order_epsilon) == false)
                    continue; // this edge should be ignore
                edge_num++;
                float tmp_sim = 0;
                for (const auto &sim : i_nei.sim_vec)
                    tmp_sim += sim;
                tmp_sim = tmp_sim / i_nei.sim_vec.size();
                community_sim += tmp_sim;
            }
        }
        community_sim = community_sim / edge_num;
        edge_all.push_back(edge_num);
        sim_all.push_back(community_sim);

        int ssc_id = query_i;
        for (auto com_i : res)
        {
            community_number[com_i] = ssc_id;
            is_in_community[com_i] = false;
            community_member[com_i] = true;
        }
        community_all_num++;
        com_size.push_back(com_num);
        sort(res.begin(), res.end());
        c_member_i.insert(c_member_i.end(), res.begin(), res.end());
        all_res_com.push_back(res);
    }

    cout << "finish find all SSC, start find hub and outlier" << endl;

    Others index_others;
    index_others.compute_hub_outlier_index(*this, c_member_i, community_number);

    t1.StopAndPrint("finished time");

    cout << "community_all_number is " << community_all_num << endl;
    cout << "all vertex :" << num_query_type_ << endl;
    int outlier_number = 0, hub_number = 0;
    int core_num = 0, community_all = 0;
    for (const auto &res : all_res_com)
    {
        community_all += res.size();
        for (const auto ci : res)
        {
            if (cand_core_[ci] == true)
                core_num++;
        }
    }
    cout << "all community number: " << community_all << endl;
    cout << "all core: " << core_num << endl;
    cout << "all non-core: " << community_all - core_num << endl;

    for (int i = 0; i < num_query_type_; i++)
    {
        if (index_others.outlier[i] == true)
            outlier_number++;
        if (index_others.hub[i] == true)
            hub_number++;
    }
    cout << "all hub: " << hub_number << endl;
    cout << "all outlier: " << num_query_type_ - community_all - hub_number << endl;
    cout << "all outlier 2: " << outlier_number - hub_number << endl;
    return;

    double cc_all = 0;
    int cnt = 0;
    long core_all = 0;
    long vertex_all = 0;
    long dia_all = 0;
    for (const auto &res : all_res_com)
    {
        int res_core = 0, res_ver = 0;
        for (const auto ci : res)
        {
            res_ver++;
            if (cand_core_[ci] == false)
                continue;
            res_core++;

            if (qn_adj_List[ci].size() == 0)
            {
                qn_adj_List[ci] = path_utils.p_induced_graph(*this, ci);
            }
        }
        core_all += res_core;
        double res_cc = global_clustering_coefficient(qn_adj_List, res);
        cc_all += res_cc;
        vertex_all += res_ver;
        int dia = compute_diameter(qn_adj_List, res);
        if (dia > 0)
            dia_all += dia;
        cnt++;
    }
    cout << "cc: is " << cc_all / cnt << endl;
    cout << "dia: is " << dia_all / cnt << endl;
    cout << "average core: is " << double(core_all / cnt) << endl;
    cout << "average vertex: is " << double(vertex_all / cnt) << endl;
    cout << "community_all_number is " << community_all_num << endl;
    int tmp_all = 0;
    long tmp_edge = 0;
    double tmp_sim = 0;
    for (auto i : com_size)
        tmp_all += i;
    for (auto i : edge_all)
        tmp_edge += i;
    for (auto i : sim_all)
        tmp_sim += i;
    double aver_size = double(tmp_all) / community_all_num;
    double aver_edge = double(tmp_edge) / community_all_num;
    double aver_sim = tmp_sim / community_all_num;
    cout << "avn: vertex number per community is " << aver_size << endl;
    cout << "ave: edge number per community is " << aver_edge << endl;
    cout << "avs: average similarity is " << aver_sim << endl;

    return;

    vector<string> vertex_names(n);
    string vertex_name_file = data_dir_ + "/vertex_name.txt";
    ifstream inputFile(vertex_name_file);
    std::string line;
    int idx = 0;
    utils u1;
    int line_cnt = u1.getLineCount(vertex_name_file);
    cout << n << " " << line_cnt << (n == line_cnt) << endl;
    while (std::getline(inputFile, line))
    {
        std::istringstream iss(line);
        int num1, num2;
        std::string name;
        iss >> num1 >> num2;
        std::getline(iss >> std::ws, name);
        while (!name.empty() && (name.back() == '\n' || name.back() == '\r'))
            name.pop_back();
        vertex_names[idx] = name;
        idx++;
    }
    cout << "finish load vertex name, start output result" << endl;

    int res_com_id = 0;
    for (const auto &res : all_res_com)
    {
        res_com_id++;
        if (res.size() > 1000 || res.size() < 10)
        {
            // cout << "this community contain too many vertex, not print" << endl;
            continue;
        }
        cout << "community num :" << res_com_id << endl;
        cout << "the number of vertex in this community: " << res.size() << endl;

        set<int> res_set(res.begin(), res.end());
        set<int> res_connect_set(res.begin(), res.end());
        int edge_all_cur = 0;
        for (auto i : res)
        {
            cout << i << " : ";
            vector<int> out_nei;
            vector<float> out_sim;
            vector<pair<int, float>> other_nei;
            out_nei.reserve(h_sim[i].size());
            out_sim.reserve(h_sim[i].size());
            other_nei.reserve(h_sim[i].size());
            for (const auto &i_nei : h_sim[i])
            {
                float sum_sim = sum_float_vec(i_nei.sim_vec);
                if (judge_demoinate(i_nei.sim_vec, index_order_epsilon) == false)
                {
                    other_nei.push_back(make_pair(i_nei.neighbor_i, sum_sim));
                    continue; // this edge should be ignore
                }
                if (res_set.find(i_nei.neighbor_i) == res_set.end())
                {
                    out_nei.push_back(i_nei.neighbor_i); // not in cur community
                    out_sim.push_back(sum_sim);          // not in cur community
                }
                else
                {
                    edge_all_cur++;
                    cout << i_nei.neighbor_i << "-" << sum_sim << " | ";
                }
            }
            cout << "***";
            if (!out_nei.empty())
            {
                for (int t = 0; t < out_nei.size(); t++)
                {
                    res_connect_set.insert(out_nei[t]);
                    cout << out_nei[t] << "-" << out_sim[t] << " | ";
                }
            }
            if (!other_nei.empty())
            {
                for (auto pair : other_nei)
                {
                    cout << pair.first << "-" << pair.second << " | ";
                    res_connect_set.insert(pair.first);
                }
            }
            cout << ": " << vertex_names[i + query_type_offset_] << endl;
        }
        cout << "current avg degree is :" << double(edge_all_cur) / res.size() << "  ";
        if (double(edge_all_cur) / res.size() > 0.4 * res.size())
            cout << " this community has high avg degree :" << res.size() << " id :" << res_com_id << endl;
        else
            cout << endl;
        for (auto i : res_connect_set)
            cout << i << " ";
        cout << endl;
    }
    cout << "finish output" << endl;

    return;

    // cout << "start save community number" << endl;

    // // string community_num_path = data_index_dir_ + "/community_num.txt";
    // string community_num_path = data_index_dir_ + "/" + to_string(p_query_type) + "/community_num" + to_string(p_mu) + ".txt";
    // ofstream community_num_file = open_file_ofstream(community_num_path);
    // for (int i = 0; i < num_query_type_; i++)
    // {
    //     if (community_number[i] == 0)
    //         continue;
    //     community_num_file << i << " " << community_number[i] << endl;
    // }
    // community_num_file.close();
}

bool HinGraph::index_judge_core(int i, int k)
{
    for (const auto &t_v : node_k_thres[i].thres_vecs)
    {
        if (judge_demoinate(t_v, index_order_epsilon))
        {
            // is_in_community[i] = true;
            return true;
        }
    }
    // is_in_community[i] = false;
    return false;
}

void HinGraph::search_k_strata(int i)
{
    if (ks_visit[i] == true)
        return;
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
            // int nei_dis = distance_[nei_type];
            // if (nei_dis == -1)
            // { // this type neighbor is unconsidered. continue
            //     continue;
            // }

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
    if (cand_sn_list[i].size() != 0)
        return;
    // 1. target type is considered and located in k-strata
    Vertex_neighbor &dn_i = dn_adj_List[i];

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
    vector<bool> v_vertex_cand(n, false);
    for (auto v_i : gen_type_dn)
    {
        v_vertex_cand[v_i] = true;
        queue<pair<int, int>> gen_q;
        gen_q.push(make_pair(v_i, 0));
        while (!gen_q.empty())
        {
            int vertex_i = gen_q.front().first, cur_step = gen_q.front().second;
            int cur_type = get_vertex_type(vertex_i);
            int cur_dis = distance_[cur_type];
            gen_q.pop();
            if (cur_step >= p_d)
                continue;
            int nei_edge_start = vertex_offset_[vertex_i];
            int nei_end = ((vertex_i + 1) == n) ? m : vertex_offset_[vertex_i + 1];
            for (int j = nei_edge_start; j < nei_end; j++)
            {
                int nei_type = edges_[j].v_type, nei_id = edges_[j].v_id;
                if (nei_type == p_query_type)
                {
                    if (v_vertex_cand[nei_id] == true)
                        continue;
                    v_vertex_cand[nei_id] = true;

                    if (mode_query == 6)
                    {
                        if (nei_id < max_scale_id)
                            cand_sn_list[i].push_back(nei_id);
                    }
                    else
                        cand_sn_list[i].push_back(nei_id);
                }
                else
                {
                    if (v_vertex_cand[nei_id] == true)
                        continue;
                    v_vertex_cand[nei_id] = true;

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

void HinGraph::check_sn(int u, queue<int> &cs_queue, queue<int> &delete_q)
{
    // 1. check core
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
        is_in_community[query_i] = false;
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

// index construct
void HinGraph::search_d_neighbor()
{
    cout << "search_all_d_neighborhood" << endl;
    cout << "all query type vertices number is " << num_query_type_ << endl;

    map<int, int> type_degree_max;
    map<int, int> type_degree_min;
    for (const auto &pair : type_epsilon)
    {
        int type = pair.first;
        type_degree_max[type] = 0;
        type_degree_min[type] = 10;
    }

    for (int i = 0; i < num_query_type_; i++)
    {
        int query_vertex_id = i + query_type_offset_;
        queue<pair<int, int>> bfs_path;
        unordered_set<int> v_vertex;
        bfs_path.push(make_pair(query_vertex_id, 0));
        v_vertex.insert(query_vertex_id);
        while (!bfs_path.empty())
        {
            int curVertex_id = bfs_path.front().first, cur_step = bfs_path.front().second;
            int cur_type = get_vertex_type(curVertex_id);
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
                // if (nei_dis == -1)
                // { // this type neighbor is unconsidered. continue
                //     continue;
                // }

                v_vertex.insert(nei_id);
                bfs_path.push(make_pair(nei_id, cur_step + 1));
                if (nei_type != p_query_type)
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
        if (i % (num_query_type_ / 10) == 0)
            cout << i << endl;
    }
    for (const auto &pair : type_epsilon)
    {
        int type = pair.first;
        cout << type << " max " << type_degree_max[type] << " min " << type_degree_min[type] << endl;
    }
    cout << endl;
}

void HinGraph::save_d_n()
{
    if (option_ == "-fidxscal")
        return;
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

void HinGraph::multi_func(int i)
{
    bool skip = true;
    for (const auto type_j : index_type_order)
    {
        vector<int> &dn_i_type_i = dn_adj_List[i].d_neighbor_[type_j];
        if (dn_i_type_i.size() > 0)
        {
            skip = false;
            break;
        }
    }
    if (skip == true)
    {
        if (i % (num_query_type_ / 100) == 0)
            cout << i << endl;
        return;
    }
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
    if (i % (num_query_type_ / 100) == 0)
        cout << i << endl;
}

void HinGraph::compute_all_similarity()
{
    cout << "compute all similarity";
    h_sim.resize(num_query_type_);

    if (num_query_type_ >= 10000)
    {
        const int num_threads = std::thread::hardware_concurrency();
        std::vector<std::thread> threads;
        for (int i = 0; i < num_query_type_; i += num_threads)
        {
            for (int j = 0; j < num_threads && (i + j) < num_query_type_; ++j)
            {
                auto thread_func = std::bind(&HinGraph::multi_func, this, i + j);
                threads.emplace_back(thread_func);
            }
            for (auto &thread : threads)
            {
                if (thread.joinable())
                {
                    thread.join();
                }
            }
            threads.clear();
        }
        return;
    }

    for (int i = 0; i < num_query_type_; i++)
    {
        bool skip = true;
        for (const auto type_j : index_type_order)
        {
            vector<int> &dn_i_type_i = dn_adj_List[i].d_neighbor_[type_j];
            if (dn_i_type_i.size() > 0)
            {
                skip = false;
                break;
            }
        }
        if (skip == true)
            continue;

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
            cout << i << endl;
    }
}

void HinGraph::save_all_similarity()
{
    cout << "start save all similarity graph" << endl;
    string type_index_path;
    if (option_ == "-fidxscal")
        type_index_path = data_index_dir_ + "/" + to_string(p_query_type) + "-" + to_string(scale_);
    else
        type_index_path = data_index_dir_ + "/" + to_string(p_query_type);
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
    string type_index_path;
    if (option_ == "-fidxscal" || option_ == "-qidxscal")
        type_index_path = data_index_dir_ + "/" + to_string(p_query_type) + "-" + to_string(scale_);
    else
        type_index_path = data_index_dir_ + "/" + to_string(p_query_type);
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
{ // return vec1 dominate vec2 or not -- True: vec1 >= vec2
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
{ // remain_type is last dim, and first type is 1 in fix_type_threshold
    fix_[fix_.size() - 1] += 1;
    for (int i = fix_.size() - 1; i >= 0; i--)
    {
        if ((fix_[i] == 101) && i != 0)
        {
            fix_[i] = 0;
            fix_[i - 1] += 1;
        }
        if (i == 0 && fix_[i] == 101)
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

void HinGraph::compute_k_threshold(int start_k)
{
    cout << "start compute k threshold" << endl;
    k_homo_graph.resize(num_query_type_);
    node_k_thres.resize(num_query_type_);
    // p_mu = 1, 2 it contain itself and dominated value.
    if (start_k == -1)
    {
        for (int k = 1; k <= 2; k++)
        {
            for (int i = 0; i < num_query_type_; i++)
            {
                if (h_sim[i].size() < k)
                    continue;
                for (auto &nei_i : h_sim[i])
                {
                    if (nei_i.domin_rank == k)
                        node_k_thres[i].thres_vecs.push_back(nei_i.sim_vec);
                }
            }
            save_k_thres_vec(k);
        }
        start_k = 3;
    }

    // compute each_threshold
    Timer t1;
    t1.Start();
    for (int k = start_k; k <= k_max; k++)
    {
        cout << k << endl;
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
            bool reload_graph = fix_[fix_.size() - 1] == 0 ? true : false;
            if (reload_graph == true)
            {
                load_all_similarity();
                cout << "current fix type and threshold is " << endl;
                for (int i = 0; i < fix_type_thresold.size(); i++)
                    cout << i << "-" << fix_type_thresold[i] << " ";
                cout << endl;
            }
            cout << "|";
            // if (fix_[0] % 30 == 0)
            //     cout << endl;
            compute_connect_k_core(k, fix_type_thresold, re_type_offset);
        }
        save_k_thres_vec(k);
    }
    t1.StopAndPrint("finsih construct index: ");
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
    for (int i = 1; i <= 101; i++)
    { // from 0.01 to 1.01
        float re_type_threshold = i / 100.0;
        vector<float> thres_sim_vec(fix_type);
        float thres_ = (i - 1) / 100.0;
        thres_sim_vec.push_back(thres_);
        bool iter_next = search_and_add_threshold(k, thres_sim_vec, re_type_threshold);
        if (iter_next == false)
            break;
    }
    return;
}

void add_new_th(const vector<float> thres_sim_vec, vector<vector<float>> &old_thres_vecs)
{
    if (old_thres_vecs.size() == 0)
    { // empty threshold
        old_thres_vecs.push_back(thres_sim_vec);
        return;
    }

    vector<vector<float>> tmp;
    tmp.reserve(old_thres_vecs.size());
    for (const auto old_vec : old_thres_vecs)
    {
        if (judge_demoinate(thres_sim_vec, old_vec))
            continue;
        tmp.push_back(old_vec); // old thres vec is not dominate by new one
    }
    if (tmp.size() == old_thres_vecs.size())
    { // new thres_sim_vec maynot add to threshold
        bool add_flag = true;
        for (const auto old_vec : old_thres_vecs)
        {
            if (judge_demoinate(old_vec, thres_sim_vec))
            {
                add_flag = false;
                break;
            }
        }
        if (add_flag == true)
            tmp.push_back(thres_sim_vec);
    }
    else
    { // new thres can add to the threshold
        tmp.push_back(thres_sim_vec);
    }
    old_thres_vecs = move(tmp);
}

bool HinGraph::search_and_add_threshold(int k, const vector<float> &thres_sim_vec, float re_type_threshold)
{
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
            if (nei_j.re_type_sim >= re_type_threshold) // only retain edge great than threshold. and cureent threshold is new thres
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
        add_new_th(thres_sim_vec, node_k_thres[cur_i].thres_vecs);
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

bool judge_geq_from_start(const vector<float> &vec1, const vector<float> &vec2, int start)
{ // return vec1[start] >= vec2 or not -- True: vec1 >= vec2
    for (int i = start; i < vec1.size(); i++)
    {
        if (vec1[i] < vec2[i])
            return false;
    }
    return true;
}

bool judge_new_edge(const vector<float> sim_vec, const vector<vector<float>> corner_points)
{
    int last = corner_points.size() - 1;
    if ((sim_vec[1] > corner_points[0][1]) || (sim_vec[0] > corner_points[last][0]))
        return true;
    for (auto &cp : corner_points)
    {
        if ((cp[1] < sim_vec[1]) && (cp[0] < sim_vec[0]))
            return true;
    }
    return false;
}

int judge_new_edge_dim1(const vector<float> sim_vec, const vector<vector<float>> corner_points)
{
    int last = corner_points.size() - 1;
    if (sim_vec[0] == 0.0 && sim_vec[1] > corner_points[0][1])
        return corner_points[0][1] * 100;
    if (sim_vec[0] > corner_points[last][0])
        return corner_points[last][1] * 100;
    bool flag_new = false;
    int min_dim1 = 100;
    for (auto &cp : corner_points)
    {
        if ((cp[1] < sim_vec[1]) && (cp[0] < sim_vec[0]))
        {
            flag_new = true;
            min_dim1 = 100 * cp[1];
        }
    }
    if (flag_new == false)
        return -1;
    return min_dim1;
}

// improve index construct
void HinGraph::improve_k_thres(int start_k)
{
    cout << "start improve compute k threshold" << endl;
    node_k_thres.resize(num_query_type_);
    // p_mu = 1, 2 it contain itself and use the dominated value.
    if (start_k == -1)
    {
        for (int k = 1; k <= 2; k++)
        {
            for (int i = 0; i < num_query_type_; i++)
            {
                if (h_sim[i].size() < k)
                    continue;
                for (auto &nei_i : h_sim[i])
                {
                    if (nei_i.domin_rank == k)
                        node_k_thres[i].thres_vecs.push_back(nei_i.sim_vec);
                }
            }
            save_k_thres_vec(k);
        }
        start_k = 3;
    }

    k_homo_graph.resize(num_query_type_, vector<k_homo_adj_node>());
    for (int i = 0; i < num_query_type_; i++)
    {
        k_homo_graph[i].reserve(h_sim[i].size());
        node_k_thres[i].used_neighbor.resize(h_sim[i].size(), false);
    }

    all_community_num = num_query_type_ / start_k;

    Timer t1;
    t1.Start();
    int k = start_k;
    cout << k << endl;
    vector<float> cons(index_type_order.size(), 0.0);
    if (cons.size() == 1)
        constraint_one_dim(k, cons, 0, 0);
    else if (cons.size() == 2)
        skyline2D(k, cons);
    else if (cons.size() == 3)
        skyline3D(k, cons);
    else
        skyline_highD(k, cons, cons.size() - 1);

    t1.StopAndPrint("finsih imporve construct index: ");
    for (int i = 0; i < num_query_type_; i++)
        sortVectorsByDimensions(node_k_thres[i].thres_vecs);
    save_k_thres_vec(k);
}

void HinGraph::skyline_highD(int k, vector<float> &cons, int cur_d)
{ // consider the current dimension under the
    if (cur_d == 2)
    {
        skyline3D(k, cons);
        return;
    }
    float d_max;
    d_max = constraint_one_dim(k, cons, cur_d, 0);
    for (int i = int(d_max * 100); i >= 0; i--)
    {
        cons[2] = i / 100.0;
        cout << cur_d << " dimension cons: " << cons[2] << endl;
        skyline_highD(k, cons, cur_d - 1);
    }
}

void HinGraph::skyline3D(int k, vector<float> &cons)
{

    float d_max;
    d_max = constraint_one_dim(k, cons, 2, 0);
    cout << d_max << endl;
    for (int i = int(d_max * 100); i >= 0; i--)
    {
        cons[2] = i / 100.0;
        cout << cons[2] << " ";
        if (i % 10 == 0)
            cout << endl;
        skyline2D(k, cons);
    }
}

bool SelectConsVertex(const vector<vector<float>> &k_threshold, const vector<float> &cons, int start_type)
{
    int end_type = cons.size();
    for (const auto &thres : k_threshold)
    {
        bool flag = true;
        for (int i = start_type; i < end_type; i++)
        {
            if (thres[i] < cons[i])
                flag = false;
        }
        if (flag == true)
            return true;
    }
    return false;
}

void HinGraph::skyline2D(int k, const vector<float> cons)
{
    bool only_2d = (cons.size() == 2);
    // 1. check any NEW_EDGE under cons(fix node). if not continue current the d1 dimension.
    bitset<SP_SIZE + 1> compute_d1_thres;
    vector<bool> fix_vertex(num_query_type_, false); // update vertex under current cons

    //  1.1 select the node under cons
    if (only_2d == false)
    {
        for (int i = 0; i < num_query_type_; i++)
            if (SelectConsVertex(node_k_thres[i].thres_vecs, cons, 2) == true)
                last_vertex[i] = true;
            else
                last_vertex[i] = false;
        // 1.2 compute the fix node and corresponding d1 threshold
        for (int i = 0; i < num_query_type_; i++)
        {
            if (last_vertex[i] == false)
                continue;
            node_k_thres[i].fix_thres_dim1.reset();
            node_k_thres[i].unsat_dim1.reset();
            bool fix_flag_node = false;
            k_homo_graph[i].clear();

            for (int j = 0; j < h_sim[i].size(); j++)
            {
                const auto &nei_i = h_sim[i][j];
                if (last_vertex[nei_i.neighbor_i] == false)
                    continue;
                if (judge_geq_from_start(nei_i.sim_vec, cons, 2) == false)
                    continue;
                k_homo_graph[i].push_back(k_homo_adj_node(nei_i.neighbor_i, 0));

                int new_fix_thres = nei_i.sim_vec[1] * 100;
                node_k_thres[i].unsat_dim1.set(new_fix_thres); // current edge not met
                int start_dim1 = judge_new_edge_dim1(nei_i.sim_vec, node_k_thres[i].corner_points);
                if (start_dim1 != -1)
                {
                    if (node_k_thres[i].used_neighbor[j] == true)
                        continue;
                    node_k_thres[i].used_neighbor[j] = true;
                    fix_flag_node = true;
                    for (int dim1 = start_dim1; dim1 <= new_fix_thres; dim1++)
                    {
                        node_k_thres[i].fix_thres_dim1.set(dim1); // maybe a new fix_thres
                    }
                }
            }
            if (fix_flag_node == true)
                fix_vertex[i] = true;
        }
    }
    else
    {
        constraint_one_dim(k, cons, 1, 0);
        constraint_one_dim(k, cons, 0, 0);
        for (int i = 0; i < num_query_type_; i++)
        {
            if (node_k_thres[i].thres_vecs.size() == 0)
                last_vertex[i] = false;
            else
                last_vertex[i] = true;
        }

        for (int i = 0; i < num_query_type_; i++)
        {
            if (last_vertex[i] == false)
                continue;
            fix_vertex[i] = true;
            node_k_thres[i].fix_thres_dim1.reset();
            node_k_thres[i].unsat_dim1.reset();
            k_homo_graph[i].clear();
            for (int dim1 = 0; dim1 < SP_SIZE + 1; dim1++)
                node_k_thres[i].fix_thres_dim1.set(dim1);
            for (auto &nei_i : h_sim[i])
            {
                if (last_vertex[nei_i.neighbor_i] == false)
                    continue;
                int new_fix_thres = nei_i.sim_vec[1] * 100;
                node_k_thres[i].unsat_dim1.set(new_fix_thres);
                k_homo_graph[i].push_back(k_homo_adj_node(nei_i.neighbor_i, 0));
            }
        }
    }

    // 1.3 prune the community without fix vertex
    int community_max = num_query_type_ / k, community_num = 1;
    vector<int> visit(num_query_type_, 0);
    vector<bool> community_del(community_max, false);
    vector<bitset<SP_SIZE + 1>> community_dim1(community_max, bitset<SP_SIZE + 1>());
    for (int i = 0; i < num_query_type_; i++)
    {
        if (visit[i] != 0)
        {
            if (community_del[visit[i]] == true)
                last_vertex[i] = false;
            else
            {
                if (fix_vertex[i] == true)
                { // current vertex is fix node, update it compute dim1
                    node_k_thres[i].fix_thres_dim1 = node_k_thres[i].fix_thres_dim1 & community_dim1[visit[i]];
                    compute_d1_thres = compute_d1_thres | node_k_thres[i].fix_thres_dim1;
                }
            }
            continue;
        }
        if (last_vertex[i] == false)
            continue;
        bitset<SP_SIZE + 1> tmp;
        bool flag_del = bfs_community(i, visit, fix_vertex, community_num, tmp);
        if (flag_del == true)
        { // current community not contain fix node. delete it
            community_del[community_num] = true;
            last_vertex[i] = false;
        }
        else
        {
            community_dim1[community_num] = move(tmp);
            if (fix_vertex[i] == true)
            { // current vertex is fix node, update it compute dim1
                node_k_thres[i].fix_thres_dim1 = node_k_thres[i].fix_thres_dim1 & community_dim1[visit[i]];
                compute_d1_thres = compute_d1_thres | node_k_thres[i].fix_thres_dim1;
            }
        }
        community_num++;
    }

    float d1_max;
    d1_max = constraint_one_dim(k, cons, 1, 0);
    vector<float> cons_1(cons);
    for (int i = int(d1_max * 100); i >= 0; i--)
    {
        if (compute_d1_thres[i] == false) // current threshold not contain new edge
            continue;
        cons_1[1] = i / 100.0;
        constraint_one_dim(k, cons_1, 0, i);
    }
}

bool HinGraph::bfs_community(int start_i, vector<int> &visit, const vector<bool> fix_vertex,
                             int community_num, bitset<SP_SIZE + 1> &c_dim1)
{
    bool flag_del = true;
    queue<int> bfs_q;
    bfs_q.push(start_i);
    visit[start_i] = community_num;
    while (!bfs_q.empty())
    {
        int cur_i = bfs_q.front();
        bfs_q.pop();
        c_dim1 = c_dim1 | node_k_thres[cur_i].unsat_dim1;
        if (flag_del == true && fix_vertex[cur_i] == true) // contain fix vertex. cannot delete
            flag_del = false;

        for (const auto &nei_i : k_homo_graph[cur_i])
        {
            if (visit[nei_i.neighbor_i] == 0)
            {
                visit[nei_i.neighbor_i] = community_num;
                bfs_q.push(nei_i.neighbor_i);
            }
        }
    }
    return flag_del;
}

float HinGraph::constraint_one_dim(int k, const vector<float> cons, int type_i, int fix_vertex_thres_)
{
    // 1. select the k-core vertex under cons
    vector<int> coreNum(num_query_type_, 0);         // core Number
    vector<bool> fix_vertex(num_query_type_, false); // update vertex under current cons
    bool allzero_ = std::all_of(cons.begin(), cons.end(), [](float val)
                                { return val == 0.0f; });
    // 1.1 select vertex satisfy cons
    if (allzero_ == false)
    {
        for (int i = 0; i < num_query_type_; i++)
        {
            if (last_vertex[i] == false)
            { // not a community vertex -- judge by the skyline2D
                community_vertex[i] = false;
                continue;
            }
            if (SelectConsVertex(node_k_thres[i].thres_vecs, cons, type_i + 1) == true)
                community_vertex[i] = true;
            else
                community_vertex[i] = false;
        }
    }
    else
    {
        vector<float> initial_corner_point(2, 0.0);
        for (int i = 0; i < num_query_type_; i++)
        {
            community_vertex[i] = true;
            node_k_thres[i].corner_points.push_back(initial_corner_point);
        }
    }

    bool flag_fix = (allzero_ == false) && (type_i == 0);

    queue<int> delete_q;
    bool contain_fix_vertex = false;
    for (int i = 0; i < num_query_type_; i++)
    {
        if (community_vertex[i] == false)
            continue;
        k_homo_graph[i].clear();
        for (const auto &nei_i : h_sim[i])
        {
            if (community_vertex[nei_i.neighbor_i] == false)
                continue;
            if (judge_demoinate(nei_i.sim_vec, cons) == false)
                continue;
            float sim_type_i = nei_i.sim_vec[type_i];
            k_homo_graph[i].push_back(k_homo_adj_node(nei_i.neighbor_i, sim_type_i));
        }
        if (flag_fix && node_k_thres[i].fix_thres_dim1[fix_vertex_thres_] == true)
        {
            fix_vertex[i] = true; // this veretex may be the fix vertex
            contain_fix_vertex = true;
        }
        coreNum[i] = k_homo_graph[i].size();
        if (coreNum[i] < k)
        {
            delete_q.push(i);
        }
    }
    if (flag_fix && contain_fix_vertex == false)
        return 0.0;

    // 1.2 core-decomposition get k-core vertex
    while (!delete_q.empty())
    {
        int cur_i = delete_q.front();
        delete_q.pop();
        community_vertex[cur_i] = false; // label not in current search
        if (flag_fix)
            fix_vertex[cur_i] = false;
        coreNum[cur_i] = 0;
        for (const auto &nei_cur_i : k_homo_graph[cur_i])
        {
            if (coreNum[nei_cur_i.neighbor_i] >= k)
            {
                coreNum[nei_cur_i.neighbor_i]--;
                if (coreNum[nei_cur_i.neighbor_i] < k)
                    delete_q.push(nei_cur_i.neighbor_i);
            }
        }
    }

    // 1.3 get connect k-core
    bitset<SP_SIZE + 2> compute_one_dim_thres;
    for (int i = 0; i < num_query_type_; i++)
    {
        if (community_vertex[i] == false)
            continue;
        vector<k_homo_adj_node> tmp_adj;
        tmp_adj.reserve(coreNum[i]);
        node_k_thres[i].one_dim.reset();
        for (const auto &nei : k_homo_graph[i])
        {
            if (community_vertex[nei.neighbor_i] == true)
            {
                tmp_adj.push_back(nei);
                int one_dim = nei.re_type_sim * 100 + 1;
                node_k_thres[i].one_dim.set(one_dim);
            }
        }
        if (tmp_adj.size() != coreNum[i])
        {
            cout << i << " error " << endl;
            cout << tmp_adj.size() << " " << coreNum[i] << endl;
            exit(-1);
        }
        compute_one_dim_thres = compute_one_dim_thres | node_k_thres[i].one_dim;
        k_homo_graph[i] = move(tmp_adj);
    }

    // ---- prune the subgraph no need to compute.
    vector<int> visit(num_query_type_, 0);
    if (flag_fix)
    {
        int community_max = num_query_type_ / k;
        vector<bool> community_del(community_max, false);
        int community_num = 1;
        compute_one_dim_thres.reset();
        for (int i = 0; i < num_query_type_; i++)
        {
            if (visit[i] != 0)
            {
                if (community_del[visit[i]] == true)
                { // current community not contain fix node. delete it
                    community_vertex[i] = false;
                }
                compute_one_dim_thres = compute_one_dim_thres | node_k_thres[i].one_dim;
                continue;
            }
            if (community_vertex[i] == false)
                continue;
            bitset<SP_SIZE + 1> tmp;
            bool flag_del = bfs_community(i, visit, fix_vertex, community_num, tmp);
            if (flag_del == true)
            { // current community not contain fix node. delete it
                community_del[community_num] = true;
                community_vertex[i] = false;
            }
            compute_one_dim_thres = compute_one_dim_thres | node_k_thres[i].one_dim;

            community_num++;
        }
        all_community_num = community_num;
    }

    // 2. compute the maximal threshold at TYPE_I dimension
    for (int i = 1; i <= 101; i++)
    {                                          // from 0.01 to 1.01
        if (compute_one_dim_thres[i] == false) // graph not contain this thres
            continue;
        float re_type_threshold = i / 100.0;
        vector<float> thres(cons);
        thres[type_i] = (i - 1) / 100.0; // add the last point to the threshold
        bool iter_next = compute_one_dim_max(k, re_type_threshold, thres,
                                             type_i, visit, fix_vertex);
        if (iter_next == false) // current max thres
            return (i - 1) / 100.0;
    }
    return (101 - 1) / 100.0;
}

bool compareCorners(const vector<float> &a, const vector<float> &b)
{
    if (a[1] > b[1]) // Compare based on the dim1, in descending order
        return true;
    else if (a[1] < b[1])
        return false;

    return a[0] > b[0]; // compare based on the dim0, in descending order
}

void update_concer_point(const vector<float> &cons, k_threshold &thres_corner, int type_i, bool fix_vertex)
// void update_concer_point(const vector<float> &cons, k_threshold &thres_corner, int type_i)
{
    // 1. update threshold
    if (thres_corner.thres_vecs.size() == 0)
        thres_corner.thres_vecs.push_back(cons);
    else
        add_new_th(cons, thres_corner.thres_vecs);

    // 2. update corner point;
    if (type_i != 0) // only update the corner point when dim0
        return;
    // if (fix_vertex)
    // {
    //     int dim_0 = cons[0] * 100;
    //     if (dim_0 >= thres_corner.max_dim0)
    //         thres_corner.fix_thres_dim1.reset();
    // }
    vector<float> new_thres(cons.begin(), cons.begin() + 2);

    vector<vector<float>> tmp_corners;
    tmp_corners.reserve(2 * thres_corner.corner_points.size());
    for (auto &&old_corner : thres_corner.corner_points)
    {
        if (old_corner[1] < new_thres[1] && old_corner[0] < new_thres[0])
        { // this corner is be dominated by the new theshold
            vector<float> add_corner1(old_corner);
            vector<float> add_corner0(old_corner);
            add_corner0[0] = new_thres[0];
            add_corner1[1] = new_thres[1];
            tmp_corners.push_back(add_corner1);
            tmp_corners.push_back(add_corner0);
        }
        else
        {
            tmp_corners.push_back(old_corner);
        }
    }
    sort(tmp_corners.begin(), tmp_corners.end(), compareCorners);
    vector<vector<float>> tmp_2;
    tmp_2.reserve(tmp_corners.size());
    for (int i = 0; i < tmp_corners.size(); i++)
    {
        bool add_flag = true;
        for (int j = i + 1; j < tmp_corners.size(); j++)
        {
            if (judge_demoinate(tmp_corners[i], tmp_corners[j]) == true)
            {
                add_flag = false;
                break;
            }
        }
        if (add_flag == true)
            tmp_2.push_back(tmp_corners[i]);
    }

    thres_corner.corner_points = move(tmp_2);
    return;
}

bool HinGraph::compute_one_dim_max(int k, float re_type_threshold, const vector<float> cons,
                                   int type_i, const vector<int> visit, vector<bool> fix_vertex)
{
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
            if (nei_j.re_type_sim >= re_type_threshold) // only retain edge great than threshold. and cureent threshold is new thres
                tmp_adj.push_back(nei_j);
        }
        k_homo_graph[i] = move(tmp_adj);
        if (k_homo_graph[i].size() < k)
        { // current vertex cannot be a core at current threshold
            delete_q.push(i);
        }
        coreNum[i] = k_homo_graph[i].size();
    }

    // 2. core decomposition
    vector<bool> community_del(all_community_num, false);
    bool fix_ = (type_i == 0);
    while (!delete_q.empty())
    {
        int cur_i = delete_q.front();
        delete_q.pop();
        community_vertex[cur_i] = false;
        if (fix_ && fix_vertex[cur_i] == true)
        { // current community may delete
            community_del[visit[cur_i]] = true;
        }
        // TODO fixnode update a new thres need to delete the fix thres.
        update_concer_point(cons, node_k_thres[cur_i], type_i, fix_vertex[cur_i] == true);
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
        if (type_i == 0 && fix_vertex[i] == true)
        { // current community should remain
            community_del[visit[i]] = false;
        }

        vector<k_homo_adj_node> tmp_adj;
        tmp_adj.reserve(coreNum[i]);
        for (auto &nei : k_homo_graph[i])
        {
            if (community_vertex[nei.neighbor_i] == true)
                tmp_adj.push_back(nei);
        }
        if (tmp_adj.size() != coreNum[i])
        {
            cout << re_type_threshold << " error " << i << endl;
            cout << tmp_adj.size() << " " << coreNum[i] << endl;
            for (const auto sim : cons)
                cout << sim << " ";
            exit(-1);
        }
        k_homo_graph[i] = move(tmp_adj);
        next_iteration = true;
    }

    bool del_c_flag = false;
    for (int i = 0; i < community_del.size(); i++)
    {
        if (community_del[i] == true)
        {
            del_c_flag = true;
            break;
        }
    }

    if (fix_ && del_c_flag)
    {
        next_iteration = false;
        for (int i = 0; i < num_query_type_; i++)
        {
            if (community_vertex[i] == false)
                continue;
            if (community_del[visit[i]] == true)
            {
                community_vertex[i] = false;
            }
            else
            {
                next_iteration = true;
            }
        }
    }
    return next_iteration;
}

void HinGraph::build_tree_index(int start_k)
{
    Timer t1;
    t1.Start();
    cout << "start build MuTree index" << endl;
    for (int j = 0; j < index_type_order.size(); j++)
    {
        index_tree[j].setMTNodes(num_query_type_);
        for (int cons = 0; cons <= 100; cons++)
        {
            build_one_type_tree(start_k, j, cons);
        }
        string mutree_index_path = data_index_dir_ + "/" + to_string(p_query_type) + "/mutree";
        index_tree[j].save_tree(mutree_index_path, p_mu, j);
        cout << "finish build MuTree for type " << index_type_order[j] << endl;
    }
    t1.StopAndPrint("finsih build MuTree index: ");
}

void HinGraph::build_one_type_tree(int k, int type_i, int cons)
{
    // 1. select the k-core vertex under cons
    vector<int> coreNum(num_query_type_, 0); // core Number
    if (cons == 0)
    {
        for (int i = 0; i < num_query_type_; i++)
        {
            k_homo_graph[i].reserve(h_sim[i].size());
            community_vertex[i] = true;
        }
    }
    double type_i_thres = cons / 100.0;
    queue<int> delete_q;
    for (int i = 0; i < num_query_type_; i++)
    {
        if (community_vertex[i] == false)
            continue;
        if (cons == 0)
        {
            k_homo_graph[i].reserve(h_sim[i].size());
            for (const auto &nei_i : h_sim[i])
                k_homo_graph[i].push_back(k_homo_adj_node(nei_i.neighbor_i, nei_i.sim_vec[type_i]));
        }
        else
        {
            vector<k_homo_adj_node> tmp_adj;
            tmp_adj.reserve(h_sim[i].size());
            for (auto nei_j : k_homo_graph[i])
            {
                if (nei_j.re_type_sim >= type_i_thres) // only retain edge great than threshold. and cureent threshold is new thres
                    tmp_adj.push_back(nei_j);
            }
            k_homo_graph[i] = move(tmp_adj);
        }
        coreNum[i] = k_homo_graph[i].size();
        if (coreNum[i] < k)
            delete_q.push(i);
    }

    // 1.2 core-decomposition get k-core vertex
    while (!delete_q.empty())
    {
        int cur_i = delete_q.front();
        delete_q.pop();
        community_vertex[cur_i] = false; // label not in current search
        coreNum[cur_i] = 0;
        for (const auto &nei_cur_i : k_homo_graph[cur_i])
        {
            if (coreNum[nei_cur_i.neighbor_i] >= k)
            {
                coreNum[nei_cur_i.neighbor_i]--;
                if (coreNum[nei_cur_i.neighbor_i] < k)
                    delete_q.push(nei_cur_i.neighbor_i);
            }
        }
    }

    // 1.3 get connect k-core
    for (int i = 0; i < num_query_type_; i++)
    {
        if (community_vertex[i] == false)
            continue;
        vector<k_homo_adj_node> tmp_adj;
        tmp_adj.reserve(coreNum[i]);
        for (const auto &nei : k_homo_graph[i])
        {
            if (community_vertex[nei.neighbor_i] == true)
                tmp_adj.push_back(nei);
        }
        if (tmp_adj.size() != coreNum[i])
        {
            cout << i << " error " << endl;
            cout << tmp_adj.size() << " " << coreNum[i] << endl;
            exit(-1);
        }
        k_homo_graph[i] = move(tmp_adj);
    }

    bool first_add = index_tree[type_i].MTNodes.size() == 1;
    // search subgraph and add to MuTree
    vector<bool> visit(num_query_type_, false);
    for (int i = 0; i < num_query_type_; i++)
    {
        if (community_vertex[i] == false)
            continue;
        if (visit[i] == true)
            continue;
        vector<int> cluster;
        cluster.reserve(initial_vector_size);
        queue<int> bfs_q;

        cluster.push_back(i);
        bfs_q.push(i);
        visit[i] = true;
        while (!bfs_q.empty())
        {
            int cur_i = bfs_q.front();
            bfs_q.pop();
            for (const auto &nei_i : k_homo_graph[cur_i])
            {
                if (visit[nei_i.neighbor_i] == false)
                {
                    cluster.push_back(nei_i.neighbor_i);
                    visit[nei_i.neighbor_i] = true;
                    bfs_q.push(nei_i.neighbor_i);
                }
            }
        }

        if (first_add)
        {
            index_tree[type_i].AddChildNode(0, type_i_thres, cluster, i);
        }
        else
        {
            index_tree[type_i].FindandUpdateNode(i, type_i_thres, cluster);
        }
    }
}

void HinGraph::intersection_neisim(vector<Nei_similarity> &nei_sim, const vector<Query_nei_similarity> &vec2)
{
    vector<Nei_similarity> intersection;
    intersection.reserve(nei_sim.size());
    int i = 0;
    int j = 0;

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
        // return 0;
        // if (n_types > 100)
        //     return 0;
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
        mergedVector.insert(mergedVector.end(), induced_qn.begin(), induced_qn.end());
    }
    sort(mergedVector.begin(), mergedVector.end());
    mVector = move(mergedVector);
}

void HinGraph::save_k_thres_vec(int k)
{
    cout << "start save threshold vector " << k << endl;
    string type_index_path;
    if (option_ == "-fidxscal")
        type_index_path = data_index_dir_ + "/" + to_string(p_query_type) + "-" + to_string(scale_);
    else
        type_index_path = data_index_dir_ + "/" + to_string(p_query_type);

    if (mode_query == 101)
    {
        type_index_path += "/improved";
    }
    check_dir_path(type_index_path);
    string k_thres_path = type_index_path + "/k_thres_" + to_string(k) + ".txt";

    ofstream k_thres_file = open_file_ofstream(k_thres_path);
    k_thres_file << num_query_type_ << " " << k << " " << index_type_order.size() << endl;
    for (int i = 0; i < num_query_type_; i++)
    {
        k_thres_file << i << " " << node_k_thres[i].thres_vecs.size() << endl;
        for (const auto &thres_vec : node_k_thres[i].thres_vecs)
        {
            for (const auto &vec_i : thres_vec)
            {
                k_thres_file << vec_i << " ";
            }
            k_thres_file << endl;
        }
        node_k_thres[i].thres_vecs.clear();
    }
    cout << "finish save threshold vector " << k << endl;
}

void HinGraph::load_k_thres_vec(int k)
{
    string type_index_path;
    if (option_ == "-qidxscal")
        type_index_path = data_index_dir_ + "/" + to_string(p_query_type) + "-" + to_string(scale_);
    else
        type_index_path = data_index_dir_ + "/" + to_string(p_query_type);

    string k_thres_path = type_index_path + "/k_thres_" + to_string(k) + ".txt";

    std::ifstream sim_file(k_thres_path.c_str());
    if (sim_file.good() == false)
    {
        k_thres_path = type_index_path + "/improved/k_thres_" + to_string(k) + ".txt";
    }

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
        node_k_thres[i].thres_vecs.resize(num_i_tv);
        for (int j = 0; j < num_i_tv; j++)
        {
            node_k_thres[i].thres_vecs[j].resize(num_index_types);
            for (int vec_i = 0; vec_i < num_index_types; vec_i++)
            {
                k_thres_file >> node_k_thres[i].thres_vecs[j][vec_i];
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
        // if(dn_adj_List[a].d_neighbor_[type].size() == 0 || dn_adj_List[b].d_neighbor_[type].size() == 0)
        //     return false;
        double sim = 0.0;
        if (option_ == "-qeffcos")
            sim = calCosSim(dn_adj_List[a].d_neighbor_[type], dn_adj_List[b].d_neighbor_[type]);
        else
            sim = calJacSim(dn_adj_List[a].d_neighbor_[type], dn_adj_List[b].d_neighbor_[type]);
        if (sim < type_ep)
            return false;
    }
    return true;
}

double HinGraph::calJacSim(const vector<int> &vec1, const vector<int> &vec2)
{
    int size_a = vec1.size(), size_b = vec2.size();
    if (size_a == 0 && size_b == 0)
    {
        return 1.0;
        // return 0.0;
    }
    else if ((size_a == 0 && size_b != 0) || (size_a != 0 && size_b == 0))
    {
        return 0.0;
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
    return static_cast<double>(intersection) / static_cast<double>(unionSize);
}

double HinGraph::calCosSim(const vector<int> &vec1, const vector<int> &vec2)
{
    int size_a = vec1.size(), size_b = vec2.size();
    if (size_a == 0 && size_b == 0)
    {
        return 1.0;
    }
    else if ((size_a == 0 && size_b != 0) || (size_a != 0 && size_b == 0))
    {
        return 0.0;
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
    double denominator = sqrt(size_a * size_b);
    return static_cast<double>(intersection) / denominator;
}

vector<double> HinGraph::avgjac_cosSim(int a, int b)
{
    if (a == b)
    {
        vector<double> vec(2, 1.0);
        return vec;
    }
    vector<double> vec(2, 0.0);
    search_k_strata(a);
    search_k_strata(b);

    double avg_jac = 0.0, avg_cos = 0.0;
    int type_num = 0;
    for (const auto &pair : type_epsilon)
    {
        int type = pair.first;
        double type_ep = pair.second;
        if (type_ep == 0.0)
            continue;
        double sim_cos = calCosSim(dn_adj_List[a].d_neighbor_[type], dn_adj_List[b].d_neighbor_[type]);
        double sim_jac = calJacSim(dn_adj_List[a].d_neighbor_[type], dn_adj_List[b].d_neighbor_[type]);
        if (sim_cos == 0 || sim_jac == 0)
            continue;
        avg_jac += sim_jac;
        avg_cos += sim_cos;
        type_num++;
    }
    if (type_num == 0)
        return vec;
    vec[0] = avg_jac / type_num;
    vec[1] = avg_cos / type_num;
    return vec;
}

double HinGraph::avgStrSim(int a, int b)
{
    if (a == b)
        return 1.0;
    double avg = 0.0;
    int type_num = 0;
    for (const auto &pair : type_epsilon)
    {
        int type = pair.first;
        double type_ep = pair.second;
        if (type_ep == 0.0)
            continue;
        double sim = 0.0;
        if (option_ == "-qeffcos")
            sim = calCosSim(dn_adj_List[a].d_neighbor_[type], dn_adj_List[b].d_neighbor_[type]);
        else
            sim = calJacSim(dn_adj_List[a].d_neighbor_[type], dn_adj_List[b].d_neighbor_[type]);
        avg += sim;
        type_num++;
    }
    return static_cast<double>(avg) / static_cast<double>(type_num);
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
            sort(empty_dn_set[type].begin(), empty_dn_set[type].end());
            int prev = empty_dn_set[type][0];
            for (int i = 1; i < empty_dn_set[type].size(); i++)
            {
                if (prev == empty_dn_set[type][i])
                    cout << prev << " ";
                else
                    prev = empty_dn_set[type][i];
            }
            cout << endl;
            cout << "type " << type << " empty query vertices contains duplicate numbers. " << endl;
        }
        cout << "type " << type << " empty query vertices number is " << empty_dn_set[type].size();
        int type_start = vertex_start_map_[type];
        int type_end = (type + 1) == n_types ? n : vertex_start_map_[type + 1];
        int max_query_vertex = 0;
        for (int j = type_start; j < type_end; j++)
        {
            if (max_query_vertex < t_hop_visited[j].size())
                max_query_vertex = t_hop_visited[j].size();
        }
        cout << " max contained " << max_query_vertex << endl;
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
    Others output_utils;
    
    output_utils.output_CD_result(*this, output);

    // string output_file_path = output + "/query_type_" + to_string(p_query_type) + "-mu" + to_string(p_mu);
    // output_file_path += "-" + data_file_name + "-" + query_file_name + "-mode" + to_string(mode_query) + ".txt";
    // ofstream output_file = open_file_ofstream(output_file_path);
    // cout << "output file name and path: " << output_file_path << endl;
    // cout << "c/n vertex_id -- cluster_id\n";

    // for (int i = 0; i < num_query_type_; i++)
    // {
    //     int vertex_id = i + query_type_offset_;
    //     if (is_in_community[i])
    //     {
    //         if (cand_core_[i])
    //         {
    //             output_file << "c " << vertex_id << endl;
    //         }
    //         else
    //         {
    //             output_file << "n " << vertex_id << endl;
    //         }
    //     }
    // }
    // output_file.close();
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

int HinGraph::find_root(int u)
{
    int x = u;
    while (pa[x] != x)
        x = pa[x];

    while (pa[u] != x)
    {
        int tmp = pa[u];
        pa[u] = x;
        u = tmp;
    }

    return x;
}

void HinGraph::my_union(int u, int v)
{
    int ru = find_root(u);
    int rv = find_root(v);

    if (ru == rv)
        return;

    if (p_rank_[ru] < p_rank_[rv])
        pa[ru] = rv;
    else if (p_rank_[ru] > p_rank_[rv])
        pa[rv] = ru;
    else
    {
        pa[rv] = ru;
        ++p_rank_[ru];
    }
}