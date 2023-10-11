#include "ONIndex.h"

ONIndex::ONIndex(/* args */)
{
}

void ONIndex::Initial_ON(int num_query_type, int type_num)
{
    all_node_size = num_query_type;
    type_size = type_num;
    ONList.resize(num_query_type);
    for (int i = 0; i < num_query_type; i++)
    {
        ONList[i].neighbor_.resize(type_num);
    }
}

bool compareBySimilarity(const Query_nei_similarity &a, const Query_nei_similarity &b)
{
    return a.similarity > b.similarity;
}

void ONIndex::save_ON(string save_path, vector<int> type_j_vec)
{
    check_dir_path(save_path);
    string on_path = save_path + "/order_neighbor.txt";
    ofstream on_file = open_file_ofstream(on_path);
    cout << "Start save the Ordered Neighbor Index. Output file name and path: " << on_path << endl;

    on_file << all_node_size << " " << type_size << endl;
    for (int i = 0; i < all_node_size; i++)
    {
        on_file << i << endl;
        for (int j = 0; j < type_size; j++)
        {
            on_file << type_j_vec[j] << " " << ONList[i].neighbor_[j].size() << " ";
            for (const Query_nei_similarity &tmp : ONList[i].neighbor_[j])
            {
                on_file << tmp.neighbor_id << " " << tmp.similarity << " ";
            }
            on_file << endl;
        }
    }
}

void ONIndex::load_ON(string load_path, vector<int> type_j_vec)
{
    cout << "start load Ordered Neighbor Index " << endl;
    string on_path = load_path + "/order_neighbor.txt";
    string line;
    ifstream on_file = open_file_fstream(on_path);

    getline(on_file, line);
    istringstream iss(line);
    iss >> all_node_size >> type_size;
    ONdelist.resize(all_node_size);
    // ONList.resize(all_node_size);
    for (int i = 0; i < all_node_size; i++)
    {
        // ONList[i].neighbor_.resize(type_size);
        ONdelist[i].neighbor_id.resize(type_size);
        ONdelist[i].neighbor_sim.resize(type_size);
    }

    int tp_num = all_node_size / 10;
    for (int i = 0; i < all_node_size; i++)
    {
        int node_id, type_j, type_j_size, id;
        getline(on_file, line);
        istringstream iss2(line);
        iss2 >> node_id;

        for (int j = 0; j < type_size; j++)
        {
            getline(on_file, line);
            istringstream iss3(line);
            iss3 >> type_j >> type_j_size;
            if (type_j != type_j_vec[j])
            {
                cout << "load ON error" << type_j << "expected type j in type_j_vec" << type_j_vec[j] << endl;
                exit(-1);
            }
            // ONList[i].neighbor_[j].reserve(type_j_size);
            ONdelist[i].neighbor_id[j].reserve(type_j_size);
            ONdelist[i].neighbor_sim[j].reserve(type_j_size);
            for (int k = 0; k < type_j_size; k++)
            {
                float similarity;
                iss3 >> id >> similarity;
                // ONList[i].neighbor_[j].emplace_back(id, similarity);
                ONdelist[i].neighbor_id[j].emplace_back(id);
                ONdelist[i].neighbor_sim[j].emplace_back(similarity);
            }
        }
        if ((i + 1) % tp_num == 0)
        {
            cout << "finish load :" << i << "/" << all_node_size << endl;
            // break;
        }
    }
    cout << "Finish load Ordered Neighbor Index " << endl;
}
