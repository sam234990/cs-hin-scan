#include "MuTree.h"

MuTree::MuTree() : root_index(0)
{
    MTNodes.emplace_back(-1.0, -1);
}

void MuTree::setMTNodes(int n)
{
    MTNodes.reserve(n);
    CoreThresholdMu.reserve(n);
}

void MuTree::AddChildNode(int parent_index, double epsilon_threshold, const vector<int> &vertices, int c_id)
{
    MTNodes.emplace_back(epsilon_threshold, vertices, c_id);
    int newNodeIndex = MTNodes.size() - 1;
    MTNodes[parent_index].children.push_back(newNodeIndex);
    // leafIndices.push_back(newNodeIndex);
}

std::vector<int> calculateDifference(const std::vector<int> &vec1, const std::vector<int> &vec2)
{
    std::vector<int> difference;

    size_t i = 0;
    size_t j = 0;

    while (i < vec1.size() && j < vec2.size())
    {
        if (vec1[i] < vec2[j])
        {
            difference.push_back(vec1[i]);
            i++;
        }
        else if (vec1[i] > vec2[j])
        {
            j++;
        }
        else
        {
            i++;
            j++;
        }
    }

    while (i < vec1.size())
    {
        difference.push_back(vec1[i]);
        i++;
    }
    // Now 'difference' contains the elements in 'vec1' that are not in 'vec2'
    return difference;
}

void MuTree::FindandUpdateNode(int c_id, double epsilon, const vector<int> &vertices)
{
    // 1. find the cluster vertices
    int node_num = FindNodeNumber(c_id);
    auto &node = MTNodes[node_num];

    vector<int> difference = calculateDifference(node.vertices, vertices);
    // 2.1 cluster unchanged, only update the epsilon
    // vertices unchanged and children is empty
    if (difference.size() == 0 && node.children.empty())
    {
        node.epsilon_threshold = epsilon;
        return;
    }
    // 2.2 cluster changed, remove duplicated verices and create new node
    node.vertices = move(difference);
    AddChildNode(node_num, epsilon, vertices, c_id);
}

int MuTree::count_all_size()
{
    long all_size = 0;
    for (auto i = 0; i < MTNodes.size(); i++)
    {
        all_size += MTNodes[i].vertices.size();
    }
    return all_size;
}

void MuTree::save_tree(string save_path, int mu, int type_j)
{
    check_dir_path(save_path);
    string mu_tree_path = save_path + "/mutree_" + to_string(mu) + "_" + to_string(type_j) + ".txt";
    ofstream mu_tree_file = open_file_ofstream(mu_tree_path);
    cout << "Start save the current MuTree. output file name and path: " << mu_tree_path << endl;

    mu_tree_file << mu << " " << MTNodes.size() << endl;
    for (int i = 0; i < MTNodes.size(); i++)
    {
        const auto &Node = MTNodes[i];
        int vertices_len = Node.vertices.size(), children_len = Node.children.size();
        mu_tree_file << i << " " << Node.epsilon_threshold << " " << Node.cluster_id << " "
                     << vertices_len << " " << children_len << endl;
        for (const auto &vert : Node.vertices)
        {
            mu_tree_file << vert << " ";
        }
        mu_tree_file << -1 << endl;
        for (const auto &child : Node.children)
        {
            mu_tree_file << child << " ";
        }
        mu_tree_file << -1 << endl;
    }
    mu_tree_file << CoreThresholdMu.size() << endl;
    for (const float &value : CoreThresholdMu)
    {
        mu_tree_file << value << endl;
    }
    mu_tree_file << -223 << endl;
    mu_tree_file.close();
}

void MuTree::load_tree(string load_path, int mu, int type_j)
{
    cout << "start load MuTree, type: " << type_j << endl;
    string MuTree_path = load_path + "/mutree_" + to_string(mu) + "_" + to_string(type_j) + ".txt";
    string line;
    ifstream mutree_file = open_file_fstream(MuTree_path);
    // getline(mutree_file, line);
    int mu_file, nodes_size;
    mutree_file >> mu_file >> nodes_size;
    if (mu_file != mu)
    {
        cout << "index file has a error with construct index " << endl;
        exit(-1);
    }
    MTNodes.resize(nodes_size);

    for (int i = 0; i < MTNodes.size(); i++)
    {
        int index, c_id, vertices_len, children_len;
        double epsilon;
        mutree_file >> index >> epsilon >> c_id >> vertices_len >> children_len;
        MTNodes[i].epsilon_threshold = epsilon;
        MTNodes[i].cluster_id = c_id;

        MTNodes[i].vertices.resize(vertices_len);
        for (int j = 0; j < vertices_len; j++)
        {
            mutree_file >> MTNodes[i].vertices[j];
        }

        int vert_end_marker;
        mutree_file >> vert_end_marker; // Read the -1 marker
        if (vert_end_marker != -1)
        {
            exit(-1);
        }

        MTNodes[i].children.resize(children_len);
        for (int j = 0; j < children_len; j++)
        {
            mutree_file >> MTNodes[i].children[j];
        }

        int child_end_marker;
        mutree_file >> child_end_marker; // Read the -1 marker
        if (child_end_marker != -1)
        {
            exit(-1);
        }
    }
    int n = 0;
    mutree_file >> n;
    CoreThresholdMu.resize(n);
    for (int i = 0; i < CoreThresholdMu.size(); i++)
    {
        mutree_file >> CoreThresholdMu[i];
    }
    int ct_end_marker;
    mutree_file >> ct_end_marker; // Read the -1 marker
    if (ct_end_marker != -223)
    {
        cout << "Error: load MuTree CoreThreshold-- Invalid file format." << endl;
        exit(-1);
    }
    mutree_file.close();
}

int MuTree::query_tree(double e_threshold, vector<vector<int>> &res, vector<int> &core_id, int query_type_offset_, vector<vector<int>> &non_core_id)
{
    // initial all the core vertices cluster subnum
    vector<int> core_cnum_id(CoreThresholdMu.size(), -2);
    for (int i = 0; i < CoreThresholdMu.size(); i++)
    {
        if (CoreThresholdMu[i] < e_threshold)
        { // vertices is not contained in any cluster is -1, non-core vertex is -3(assigned later)
            core_cnum_id[i] = -1;
        }
    }

    // 1. from the root BFS the cluster start node
    vector<int> cluster_start_nodes;
    cluster_start_nodes.reserve(MTNodes.size());
    queue<int> q;
    q.push(0);
    while (!q.empty())
    {
        int current_node_index = q.front();
        q.pop();

        for (int child_index : MTNodes[current_node_index].children)
        {
            const MTNode &child_node = MTNodes[child_index];
            if (child_node.epsilon_threshold >= e_threshold)
            {
                cluster_start_nodes.push_back(child_index);
            }
            else
            {
                q.push(child_index);
            }
        }
    }

    // 2. from the cluster start node to the end.
    int cluster_num = cluster_start_nodes.size();
    vector<vector<int>> cluster_result;
    vector<int> cluster_id;
    cluster_result.resize(cluster_num);
    cluster_id.resize(cluster_num);
    for (int i = 0; i < cluster_num; i++)
    {
        int start_node_index = cluster_start_nodes[i];
        cluster_result[i].reserve(100);
        cluster_id[i] = MTNodes[start_node_index].cluster_id;

        queue<int> cluster_queue;
        cluster_queue.push(start_node_index);
        while (!cluster_queue.empty())
        {
            int current_node_index = cluster_queue.front();
            cluster_queue.pop();
            const MTNode &current_node = MTNodes[current_node_index];
            // for(const auto &x :current_node.vertices){}
            cluster_result[i].insert(cluster_result[i].end(), current_node.vertices.begin(), current_node.vertices.end());

            for (auto &child_index : current_node.children)
            {
                cluster_queue.push(child_index);
            }
        }
        sort(cluster_result[i].begin(), cluster_result[i].end());
        cluster_result[i].erase(unique(cluster_result[i].begin(), cluster_result[i].end()), cluster_result[i].end());
    }
    int core_num = 0;
    // TODO add a non-core cluster id map
    vector<vector<int>> non_core_cnum_id(CoreThresholdMu.size(), vector<int>());
    for (int i = 0; i < cluster_num; i++)
    {
        const auto &clu = cluster_result[i];
        for (const auto &x : clu)
        {
            int a = x - query_type_offset_;
            if (core_cnum_id[a] == -1)
            { // non-core in cluster is -3
                core_cnum_id[a] = NONCORE_CLUSTER_FLAG;
                non_core_cnum_id[a].push_back(i);
                continue;
            }
            else if (core_cnum_id[a] == NONCORE_CLUSTER_FLAG)
            { // non-core vertex and set FLAG before
                non_core_cnum_id[a].push_back(i);
                continue;
            }
            else
            {
                if (core_cnum_id[a] == -2)
                { // core vertex
                    core_num++;
                    core_cnum_id[a] = i;
                }
                else
                {
                    cout << "core cluster number error" << endl;
                }
            }
        }
    }
    for (int i = 0; i < CoreThresholdMu.size(); i++)
    {
        if (core_cnum_id[i] == -2)
        {
            cout << "error core: " << i << endl;
            PrintFindAllNodeNumber(i); // 1849 1947
            PrintDFSNode(FindNodeNumber(i));
        }
    }
    res = move(cluster_result);
    core_id = move(core_cnum_id);
    non_core_id = move(non_core_cnum_id);
    cout << "finish query mutree type" << endl;
    return core_num;

    // for debug
    /*
    int allvertices = 0;
    vector<pair<int, int>> result;
    for (int i = 0; i < cluster_num; i++)
    {
        const auto &i_clu = cluster_result[i];
        allvertices += i_clu.size();
        for (const auto cluster_i : i_clu)
        {
            result.push_back(make_pair(cluster_i, cluster_id[i]));
        }
    }
    cout << "clusters number is " << cluster_num << endl;
    cout << "all clusters' vertices number is " << allvertices << endl;
    sort(result.begin(), result.end());
    for (int i = 1; i < result.size(); ++i)
    {
        if (result[i] == result[i - 1])
        {
            cout << "the first duplicated result " << result[i].first << " -- " << result[i].second << endl; // 2147 -- 452
            PrintFindAllNodeNumber(result[i].first);                                                         // 1849 1947
            PrintDFSNode(FindNodeNumber(result[i].first));
            break;
        }
    }

    result.erase(unique(result.begin(), result.end()), result.end());
    cout << "all clusters' vertices number is " << result.size() << endl;
    for (const auto pair : result)
    {
        cout << pair.first << "--" << pair.second << endl;
    }
    cout << "all clusters' vertices number is " << result.size() << endl;
    */
}

void MuTree::query_tree_only_node(double e_threshold, vector<int> &vertex_res, int query_type_offset_)
{
    vector<int> core_cnum_id(CoreThresholdMu.size(), -4);
    
    // 1. from the root BFS the cluster start node
    vector<int> cluster_start_nodes;
    cluster_start_nodes.reserve(MTNodes.size());
    queue<int> q;
    q.push(0);
    while (!q.empty())
    {
        int current_node_index = q.front();
        q.pop();

        for (int child_index : MTNodes[current_node_index].children)
        {
            const MTNode &child_node = MTNodes[child_index];
            if (child_node.epsilon_threshold >= e_threshold)
            {
                cluster_start_nodes.push_back(child_index);
            }
            else
            {
                q.push(child_index);
            }
        }
    }

    // 2. from the cluster start node to the end.
    int cluster_num = cluster_start_nodes.size();
    for (int i = 0; i < cluster_num; i++)
    {
        int start_node_index = cluster_start_nodes[i];
        queue<int> cluster_queue;
        cluster_queue.push(start_node_index);
        while (!cluster_queue.empty())
        {
            int current_node_index = cluster_queue.front();
            cluster_queue.pop();
            const MTNode &current_node = MTNodes[current_node_index];
            for(const auto &vertex_id_x :current_node.vertices)
            {
                int vertex_i = vertex_id_x - query_type_offset_;
                core_cnum_id[vertex_i] = NONCORE_CLUSTER_FLAG;
            }
            for (auto &child_index : current_node.children)
            {
                cluster_queue.push(child_index);
            }
        }
    }
    for (int i = 0; i < CoreThresholdMu.size(); i++)
    {
        if (CoreThresholdMu[i] < e_threshold)
        { // vertices is not contained in any cluster is -1, non-core vertex is -3(assigned later)
            // core_cnum_id[i] = -1;
            continue;
        }
        else
        {// core vertex is assign as 0
            core_cnum_id[i] = 0;
        }
    }
    vertex_res = move(core_cnum_id);
    return;
}

int MuTree::FindNodeNumber(int c_id)
{
    for (auto i = 0; i < MTNodes.size(); i++)
    {
        const auto &Node = MTNodes[i];
        for (const auto &ver : Node.vertices)
        {
            if (c_id == ver)
            {
                return i;
            }
        }
    }
    return 0;
}

void MuTree::PrintFindAllNodeNumber(int c_id)
{
    for (auto i = 0; i < MTNodes.size(); i++)
    {
        const auto &Node = MTNodes[i];
        for (const auto &ver : Node.vertices)
        {
            if (c_id == ver)
            {
                cout << i << " ";
                break;
            }
        }
    }
    cout << endl;
}

void MuTree::PrintDFSPaths(int start_index, vector<int> &current_path)
{
    const MTNode &current_node = MTNodes[start_index];
    current_path.push_back(start_index);

    if (current_node.children.empty())
    {
        for (const int &node_index : current_path)
        {
            cout << node_index << "-" << MTNodes[node_index].epsilon_threshold << "-->";
        }
        cout << endl;
    }
    else
    {
        for (const int &child_index : current_node.children)
        {
            PrintDFSPaths(child_index, current_path);
        }
    }

    current_path.pop_back();
}

void MuTree::PrintDFSNode(int start_index)
{
    vector<int> current_path;
    PrintDFSPaths(start_index, current_path);
}