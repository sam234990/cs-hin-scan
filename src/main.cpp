#include "utils.h"
#include "HinGraph.h"

int main(int argc, char const *argv[])
{
    if (argc < 4)
    {
        cerr << "Invalid number of arguments." << endl;
        cerr << "Usage: \n"
             << argv[0] << " -f <input_path> <output_path> \n"
             << argv[0] << " -q <input> <query_file> [-o <output_path>]\n"
             << argv[0] << " -q[0-9] <input> <query_file> [-o <output_path>]"
             << endl;
        return 1;
    }

    string option = argv[1];

    if (option == "-f")
    { // process raw graph
        string input_path = argv[2];
        string output_path = argv[3];
        cout << "Option: -f" << endl;
        cout << "Input Path: " << input_path << endl;
        cout << "Output Path: " << output_path << endl;
        utils u1;
        u1.precoss_graph(input_path, output_path);
        get_memory_usage();
    }
    else if (option[1] == 'q')
    { // query
        if (argc < 4)
        {
            cerr << "Invalid number of arguments." << endl;
            cerr << "Usage: " << argv[0] << " -q[0-9 | idx] <input> <query_file> [-o <output_path>]"
                 << endl;
            return 1;
        }
        string input_path = argv[2];
        string query_file = argv[3];
        string output_path = ""; // 默认为空
        if (argc > 5)
        {
            output_path = argv[5];
        }
        cout << "Option: " << option << endl;
        cout << "Input Path: " << input_path << endl;
        cout << "Query File: " << query_file << endl;
        HinGraph g1 = HinGraph(input_path);

        g1.load_graph();
        g1.cs_hin_scan(query_file, option);

        if (!output_path.empty())
        {
            cout << "Output Path: " << output_path << endl;
            g1.output_result(output_path);
        }
        else
        {
            cout << "Output Path not provided. " << endl;
        }
        get_memory_usage();
    }
    else if (option == "-fidx" || option == "-fidxON")
    { // build index
        if (argc < 4)
        {
            cerr << "Invalid number of arguments." << endl;
            cerr << "Usage: " << argv[0] << " -fidx <input> <idx_query_file>" << endl;
            return 1;
        }
        string input_path = argv[2];
        string idx_query_file = argv[3];
        int start_k = std::stoi(argv[4]);
        cout << "Option: -fidx" << endl;
        cout << "Input Path: " << input_path << endl;
        cout << "Index Query File: " << idx_query_file << endl;
        HinGraph g1 = HinGraph(input_path);
        g1.load_graph();
        // Call the function to build the index here
        // ...
        g1.construct_index(idx_query_file, option, start_k);
        get_memory_usage();
    }
    else
    {
        cerr << "Invalid option: " << option << endl;
        return 1;
    }

    return 0;
}
