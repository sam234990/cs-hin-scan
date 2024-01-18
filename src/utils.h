#ifndef _UTILS_H_
#define _UTILS_H_

#include <ctime>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <string>
#include <cstring>
#include <vector>
#include <map>
#include <set>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#define _LINUX_

#ifdef _LINUX_
#include <sys/time.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/resource.h>
#else
#include <io.h>
#include <direct.h>
#include <windows.h>
#endif

using namespace std;

FILE *open_file(const char *file_name, const char *mode);
ifstream open_file_fstream(string &file_name);
ofstream open_file_ofstream(string &file_name);
void get_memory_usage();
void check_dir_path(string dir_path);


const int MAX_ID = 1000000000; // 10E9 1B

struct Vertex_type
{
	int v_id;
	int v_type;
	int edge_type;
};

class utils
{
	string input_dir_;
	string output_dir_;
	string dest_dir_;
	int m, n, n_types;
	int *vertex_;
	// vector<int> vertex_offset_;
	vector<int> vertex_start_map_;
	map<int, int> s_r_vertex_map_;
	vector<vector<Vertex_type>> edges_;
	// Vertex_type *edges_; // hin graph
	vector<vector<int>> hin_schema_adjacencyMatrix;
	vector<vector<int>> hin_schema_edge_cnt;

	int getLineCount(string &filename);
	int get_vertex_type(int vertex_id);
	bool static compareVertex(const Vertex_type& v1, const Vertex_type& v2);
	void make_hin_schema();
	void save_graph();
	void process_vertex();

public:
	vector<int> edge_type_;
	utils();
	~utils();
	void precoss_graph(string input, string output);
	void load_edge_type(string input);

	void read_and_print_vertex();
	void read_and_print_graph();
	void print_file(string input);
};

#endif