#include "utils.h"
#include "HinGraph.h"

int main(int argc, char const *argv[])
{
    if (argc < 4)
    {
        if (argc == 2)
        {
            string option = argv[1];
            if (option == "-help_query")
            {
                string query_schema = "Here is the schema of query file in each line:\n\n"
                                      "mu -- parameter mu"
                                      "k -- k-multi-strata \n"
                                      "qt numq -- query type & the number of query vertices\n"
                                      "vertices [-1 | -2] -- the numq query vertices  [random query(-1) | selected_file(-2)] \n"
                                      "dimension(type) value [-1 | -2 | -3 | -4] -- vector of epsilon in Jaccard measure or  other methods:\n"
                                      "\t unit value (-1)\n"
                                      "\t PathSim(-2)\n"
                                      "\t Jaccard query with PathSim metric(-3)\n"
                                      "\t PathSim query with Jaccard metric(-4)\n"
                                      "\n For other query instruction, please refer to function  load_query_file()  in HinGraph.cpp\n";

                cout << query_schema << endl;
                return 1;
            }
        }
        cerr << "Invalid number of arguments." << endl;
        cerr << "Usage: \n"
             << argv[0] << " -f <input_path> <output_path> \n"
             << argv[0] << " -q[qcdSCAN | qscal] <input> <query_file> [-o <output_path>] [scale]\n"
             << argv[0] << " -qidx[qidx] <input> <query_file> [-o <output_path>] [scale]\n"
             << argv[0] << " -fidx[fidx1 | fidxscal] <input> <output_path> [scale]\n"
             << argv[0] << " -meta <input> type\n"
             << argv[0] << " -transhomo <input> <query_file> -o <output_path> \n";
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
            cerr << "Usage: " << argv[0] << " -q[0-9 | idx] <input> <query_file> [-o <output_path>] [scale]"
                 << endl;
            return 1;
        }
        string input_path = argv[2];
        string query_file = argv[3];
        string output_path = ""; 
        int scale = 100;
        if (option == "-qidxscal" || option == "-qscal")
        {
            scale = std::stoi(argv[4]);
            cout << "scale test:" << scale << endl;
        }
        else if (argc > 5)
        {
            output_path = argv[5];
        }
        cout << "Option: " << option << endl;
        cout << "Input Path: " << input_path << endl;
        cout << "Query File: " << query_file << endl;
        HinGraph g1 = HinGraph(input_path);

        g1.load_graph();
        g1.cs_hin_scan(query_file, option, scale);

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
    else if (option == "-fidx" || option == "-fidx1" || option == "-fidxscal" || option == "-fidx1scal")
    { // build index
        if (argc < 4)
        {
            cerr << "Invalid number of arguments." << endl;
            cerr << "Usage: " << argv[0] << " -fidx <input> <idx_query_file> [scal]" << endl;
            return 1;
        }
        string input_path = argv[2];
        string idx_query_file = argv[3];
        int start_k = std::stoi(argv[4]);
        cout << "Option: " << option << endl;
        cout << "Input Path: " << input_path << endl;
        cout << "Index Query File: " << idx_query_file << endl;
        HinGraph g1 = HinGraph(input_path);
        g1.load_graph();
        // Call the function to build the index here
        // ...
        int scale = 100;
        if (option == "-fidxscal" || option == "-fidx1scal")
        {
            scale = std::stoi(argv[5]);
            cout << "scale test:" << scale << endl;
        }
        g1.construct_index(idx_query_file, option, start_k, scale);
        get_memory_usage();
    }
    else if (option == "-meta")
    {
        if (argc < 4)
        {
            cerr << "Invalid number of arguments." << endl;
            cerr << "Usage: " << argv[0] << " -meta <input> type" << endl;
            return 1;
        }
        string input_path = argv[2];
        int query_type = std::stoi(argv[3]);

        cout << "Option: " << option << endl;
        cout << "Input Path: " << input_path << endl;
        cout << "meta-path type: " << query_type << endl;

        HinGraph g1 = HinGraph(input_path);
        g1.load_graph();
        g1.find_meta(query_type);
        get_memory_usage();
    }
    else if (option == "-transhomo")
    {
        if (argc < 4)
        {
            cerr << "Invalid number of arguments." << endl;
            cerr << "Usage: " << argv[0] << " -q[0-9 | idx] <input> <query_file> [-o <output_path>] [scale]"
                 << endl;
            return 1;
        }
        string input_path = argv[2];
        string query_file = argv[3];
        string output_path = "";
        if (argc > 5)
        {
            output_path = argv[5];
        }
        else
        {
            cerr << "Invalid number of arguments." << endl;
            cerr << "Usage: " << argv[0] << " -transhomo <input> <query_file> -o <output_path> " << endl;
            return 1;
        }
        cout << "Option: " << option << endl;
        cout << "Input Path: " << input_path << endl;
        cout << "Query File: " << query_file << endl;
        HinGraph g1 = HinGraph(input_path);

        g1.load_graph();
        g1.load_query_file(query_file);
        PathSim p1 = PathSim();
        p1.trans_homo_graph(g1, g1.metapath_vecs[0], output_path);
    }
    else
    {
        cerr << "Invalid option: " << option << endl;
        return 1;
    }

    return 0;
}
