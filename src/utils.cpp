#include "utils.h"

FILE *open_file(const char *file_name, const char *mode)
{
	FILE *f = fopen(file_name, mode);
	if (f == NULL)
	{
		printf("Can not open file: %s\n", file_name);
		exit(1);
	}

	return f;
}

ifstream open_file_fstream(string &file_name)
{
	ifstream file(file_name);
	if (!file.is_open())
	{
		cerr << "Can not open file: " << file_name << endl;
		exit(1);
	}

	return file;
}
ofstream open_file_ofstream(string &file_name)
{
	ofstream output_file(file_name);
	if (!output_file.is_open())
	{
		std::cerr << "Failed to open file." << std::endl;
		exit(1);
	}
	return output_file;
}

void format_memory_usage(long memory_used)
{
	if (memory_used >= 1024 * 1024 * 1024)
	{
		double gb = static_cast<double>(memory_used) / (1024 * 1024 * 1024);
		std::cout << "Memory usage: " << std::fixed << std::setprecision(2) << gb << " GB\n";
	}
	else if (memory_used >= 1024 * 1024)
	{
		double mb = static_cast<double>(memory_used) / (1024 * 1024);
		std::cout << "Memory usage: " << std::fixed << std::setprecision(2) << mb << " MB\n";
	}
	else if (memory_used >= 1024)
	{
		double kb = static_cast<double>(memory_used) / 1024;
		std::cout << "Memory usage: " << std::fixed << std::setprecision(2) << kb << " KB\n";
	}
	else
	{
		std::cout << "Memory usage: " << memory_used << " bytes\n";
	}
}

void get_memory_usage()
{
#ifdef _WIN32
	PROCESS_MEMORY_COUNTERS_EX pmc;
	GetProcessMemoryInfo(GetCurrentProcess(), (PROCESS_MEMORY_COUNTERS *)&pmc, sizeof(pmc));
	SIZE_T memory_used = pmc.WorkingSetSize;
	format_memory_usage(memory_used);
#elif defined(_LINUX_)
	struct rusage r_usage;
	getrusage(RUSAGE_SELF, &r_usage);
	long memory_used = r_usage.ru_maxrss * 1024; // Convert from KB to bytes
	format_memory_usage(memory_used);
#endif
}

void check_dir_path(string dir_path)
{
	// make destDir directory
	if (access(dir_path.c_str(), 0) == -1)
	{
// use for LINUX
#ifdef _LINUX_
		int flag = mkdir(dir_path.c_str(), 0777);
#else
		int flag = mkdir(destdir_path_dir_.c_str());
#endif
		if (flag == 0)
			printf("Directory \"%s\" is created\n", dir_path.c_str());
		else
		{
			printf("Create directory failed\n");
			exit(1);
		}
	}
}

utils::utils()
{
	m = 0;
	n = 0;
	n_types = 0;
}

utils::~utils()
{
	// if (vertex_ != NULL)
	// 	delete[] vertex_;
	// if (edges_ != NULL)
	// 	delete[] edges_;
	// if (vertex_offset_ != NULL)
	// 	delete[] vertex_offset_;
}

void utils::precoss_graph(string input, string output)
{
	input_dir_ = input;
	output_dir_ = output;
	dest_dir_ = output + "/pro_graph";
	check_dir_path(dest_dir_);
	process_vertex();
	load_edge_type(input);

	// read_and_print_vertex();
	read_and_print_graph();
	make_hin_schema();
	save_graph();
	// delete[] vertex_;
	// delete[] edges_;
	// delete[] vertex_offset_;
	return;
}

void utils::load_edge_type(string input)
{
	string edge_file_path = input + "/edge.txt";
	int line_cnt = getLineCount(edge_file_path);
	edge_type_.resize(line_cnt);
	ifstream edge_file = open_file_fstream(edge_file_path);
	int id, type;
	for (int i = 0; i < line_cnt; i++)
	{
		edge_file >> id >> type;
		if (id != i)
		{
			cerr << "error: wrong edge file format" << endl;
			exit(-1);
		}
		edge_type_[i] = type;
	}
}

void utils::read_and_print_vertex()
{
	string vertex_file_path = input_dir_ + "/vertex.txt";

	string line, last_line;
	ifstream vertex_file = open_file_fstream(vertex_file_path);

	// save every line read from txtFile
	vertex_ = new int[MAX_ID];
	int vertex_id, type_id, cur_type_id = -1, num = 0;

	printf("Read Vertex: \n");
	while (getline(vertex_file, line))
	{
		if (line[0] < '0' || line[0] > '9')
			continue;
		sscanf(line.c_str(), "%d\t%d", &vertex_id, &type_id);
		vertex_[vertex_id] = vertex_id;
		if (type_id != cur_type_id)
		{
			cur_type_id = type_id;
			// cout << line << endl;
			vertex_start_map_.push_back(vertex_id);
		}
		// last_line = line;
	}

	for (int i : vertex_start_map_)
	{
		cout << i << " ";
	}
	cout << endl;

	n_types = vertex_start_map_.size();
	n = vertex_id + 1;
	cout << "n and n_types " << n << " " << n_types << endl;

	vertex_file.close();
	return;
}

void utils::read_and_print_graph()
{
	// count edges
	string edge_file = input_dir_ + "/edge.txt";
	int line_num = getLineCount(edge_file);
	cout << line_num << endl;
	m = line_num;

	// vertex_offset_ = new int[MAX_ID];
	edges_.resize(n);

	string graph_file_path = input_dir_ + "/graph.txt";
	string line;
	int vertex_id = -1, neighbor_vertex_id = -1, vertex_edge_num = 0;
	ifstream graph_file = open_file_fstream(graph_file_path);
	cout << "read HIN graph" << endl;
	while (getline(graph_file, line))
	{
		istringstream iss(line);
		int source_id;
		iss >> source_id;
		vertex_id = s_r_vertex_map_[source_id];

		int neighbor_source_id, edge_id;
		while (iss >> neighbor_source_id)
		{
			neighbor_vertex_id = s_r_vertex_map_[neighbor_source_id];
			Vertex_type neighbor;
			neighbor.v_id = neighbor_vertex_id;
			neighbor.v_type = get_vertex_type(neighbor_vertex_id);
			// int edge_id;
			if (!(iss >> edge_id))
			{
				break;
			}
			neighbor.edge_type = edge_type_[edge_id];
			edges_[vertex_id].push_back(neighbor);
		}
		vertex_edge_num++;
		if (vertex_edge_num % 500000 == 0)
		{
			cout << vertex_edge_num << " ";
			int cur_end = edges_[vertex_id].size() - 1;
			if (cur_end < 0)
				cout << "-1 -1" << endl;
			else
				cout << edges_[vertex_id][cur_end].v_id << " " << edges_[vertex_id][cur_end].v_type << endl;
		}
	}

	// sort each neighbor
	int start = 0, end = 0, edge_cnt = 0;

	for (int i = 0; i < n; i++)
	{
		sort(edges_[i].begin(), edges_[i].end(), compareVertex);
		// auto it = std::unique(edges_[i].begin(), edges_[i].end(),
		// 					  [](const Vertex_type &v1, const Vertex_type &v2)
		// 					  {
		// 						  return v1.v_id == v2.v_id;
		// 					  });
		// edges_[i].erase(it, edges_[i].end());
		edge_cnt += edges_[i].size();
	}
	m = edge_cnt;
	cout << m << " -- " << line_num << " -- " << int(m == line_num) << endl;
	graph_file.close();
	return;
}

void utils::print_file(string input)
{
	cout << "print the first 50 line : " << input << endl;
	ifstream graph_file = open_file_fstream(input);
	string line;
	int line_num = 0;
	cout << 1 << endl;
	cout << getLineCount(input) << endl;
	while (getline(graph_file, line))
	{
		cout << line << endl;
		line_num++;
		if (line_num > 50)
			break;
	}
}

bool utils::compareVertex(const Vertex_type &v1, const Vertex_type &v2)
{
	if (v1.v_type < v2.v_type)
	{
		return true;
	}
	else if (v1.v_type == v2.v_type)
	{
		return v1.v_id < v2.v_id;
	}
	else
	{
		return false;
	}
}

void utils::make_hin_schema()
{
	hin_schema_adjacencyMatrix.resize(n_types, std::vector<int>(n_types, 0));
	hin_schema_edge_cnt.resize(n_types, vector<int>(n_types, 0));
	for (int i = 0; i < n_types; i++)
	{
		int start_v = vertex_start_map_[i];
		int end_v = ((i + 1) == n_types) ? n : vertex_start_map_[i + 1];
		set<int> uniqueTypes;
		map<int, int> type_edge_cnt;
		for (int j = start_v; j < end_v; j++)
		{
			for (auto nei : edges_[j])
			{
				uniqueTypes.insert(nei.v_type);
				type_edge_cnt[nei.v_type]++;
			}
		}
		for (int num : uniqueTypes)
		{
			if (hin_schema_adjacencyMatrix[i][num] == 0)
			{
				hin_schema_adjacencyMatrix[i][num] = 1;
				hin_schema_adjacencyMatrix[num][i] = 1;
			}
		}
		for (auto pair : type_edge_cnt)
		{
			hin_schema_edge_cnt[i][pair.first] = pair.second;
		}
	}
	for (const auto &row : hin_schema_adjacencyMatrix)
	{
		for (int val : row)
		{
			std::cout << val << " ";
		}
		std::cout << std::endl;
	}
	cout << "edge cnt" << endl;
	for (const auto &row : hin_schema_edge_cnt)
	{
		for (int val : row)
		{
			std::cout << val << " ";
		}
		std::cout << std::endl;
	}
	return;
}

void utils::save_graph()
{
	cout << "start save graph " << endl;
	// graph info
	string graph_info_path = dest_dir_ + "/graph_info.txt";
	ofstream graph_info = open_file_ofstream(graph_info_path);
	graph_info << m << " " << n << " " << n_types << endl;
	for (int i : vertex_start_map_)
	{
		graph_info << i << " ";
	}
	graph_info << endl;
	for (int i = 0; i < n_types; i++)
	{
		for (int j = 0; j < n_types; j++)
		{
			graph_info << hin_schema_adjacencyMatrix[i][j] << " ";
		}
		graph_info << endl;
	}
	for (int i = 0; i < n_types; i++)
	{
		for (int j = 0; j < n_types; j++)
		{
			graph_info << hin_schema_edge_cnt[i][j] << " ";
		}
		graph_info << endl;
	}
	graph_info.close();

	// hin graph edges neighbor id & type
	string graph_type_path = dest_dir_ + "/graph_type.txt";
	ofstream graph_type = open_file_ofstream(graph_type_path);
	int cur_edge = 0;
	vector<int> vertex_edge_offset;
	vertex_edge_offset.resize(n);
	for (int i = 0; i < n; i++)
	{
		vertex_edge_offset[i] = cur_edge;
		const vector<Vertex_type> neighborhood_set = edges_[i];
		for (auto nei : neighborhood_set)
		{
			graph_type << nei.v_id << " " << nei.v_type << " " << nei.edge_type << endl;
			cur_edge++;
		}
	}
	graph_type.close();
	int line_num = getLineCount(graph_type_path);
	cout << line_num << " " << m << " " << int(line_num == m) << endl;

	// vertex offset
	string vertex_offset_path = dest_dir_ + "/graph_offset.txt";
	ofstream vertex_offset = open_file_ofstream(vertex_offset_path);
	for (int i = 0; i < n; i++)
	{
		vertex_offset << vertex_edge_offset[i] << endl;
	}
	vertex_offset.close();
	line_num = getLineCount(vertex_offset_path);
	cout << line_num << " " << n << " " << int(line_num == n) << endl;
	cout << "save graph finish" << endl;
}

void utils::process_vertex()
{
	string vertex_file_path = input_dir_ + "/vertex.txt";
	ifstream vertex_file = open_file_fstream(vertex_file_path);
	string line;
	map<int, vector<int>> vertex_map;
	int line_count = 0;
	while (getline(vertex_file, line))
	{
		int id, type;
		istringstream iss(line);
		iss >> id >> type;
		vertex_map[type].push_back(id);
		if (line_count++ % 500000 == 0)
			cout << line_count << endl;
	}
	vertex_file.close();
	n_types = vertex_map.size();

	// save vertex map : resorted_id source_id resorted_type
	string vertex_map_path = dest_dir_ + "/vertex_map.txt";
	ofstream vertex_map_file = open_file_ofstream(vertex_map_path);
	int cur_v_num = 0;
	int i = 0;
	for (auto it = vertex_map.begin(); it != vertex_map.end(); it++, i++)
	{
		vertex_start_map_.push_back(cur_v_num);
		cout << "type : " << i << " -- start_num : " << cur_v_num << " source type: " << it->first << endl;
		const vector<int> v_st = it->second;
		for (int source_id : v_st)
		{
			vertex_map_file << cur_v_num << " " << source_id << " " << i << endl;
			s_r_vertex_map_[source_id] = cur_v_num;
			cur_v_num++;
		}
	}
	vertex_map_file.close();
	n = cur_v_num;
}

int utils::getLineCount(string &filename)
{
	ifstream file = open_file_fstream(filename);
	int lineCount = 0;
	string line;
	while (getline(file, line))
	{ // 逐行读取文件内容
		lineCount++;
	}

	file.close(); // 关闭文件
	return lineCount;
}

int utils::get_vertex_type(int vertex_id)
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
