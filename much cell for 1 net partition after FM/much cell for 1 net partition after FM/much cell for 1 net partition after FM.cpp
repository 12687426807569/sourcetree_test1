using namespace std;
constexpr int max_enumerate_cell_num = 20;
#define TEST_CASE_NUM 1
#include <omp.h>
#include <iostream>
#include<fstream>
#include<numeric>
#include <cfloat>
#include<tuple>
//#include <map>
#include <vector>
#include <string>
#include<climits>
//#include<set>
#include<unordered_set>
#include<unordered_map>
//#include <list>
#include<cmath>
#include <algorithm>
#include <random>
#include<ctime>
#include <array>
#include<chrono>
#include <memory>
#define locked true
#define unlocked false
#define A   0
#define B   1
using cell_input_id_type = int;
using cell_program_id_type = int;
typedef int net_id_type;
using partition_name_type = int;
class IRandomIntGenerator {
public:
	virtual int Next(int min, int max) = 0;
	virtual ~IRandomIntGenerator() = default;
};
class Mt19937Generator : public IRandomIntGenerator {
	std::mt19937 rng;
public:
	Mt19937Generator() : rng(static_cast<unsigned int>(time(nullptr))) {}
	Mt19937Generator(unsigned int seed) : rng(seed) {}
	int Next(int min, int max) override {
		std::uniform_int_distribution<int> dist(min, max);
		return dist(rng);
	}
};
class RandModGenerator : public IRandomIntGenerator {
public:
	RandModGenerator() {
		srand(static_cast<unsigned int>(time(nullptr)));
	}
	RandModGenerator(unsigned int seed) {
		srand(seed);
	}
	int Next(int min, int max) override {
		return rand() % (max - min + 1) + min;
	}
};
struct myGraph {
	vector<int>cell_program_id_to_area;//V[cell_program_id]=area
	vector<cell_input_id_type>cell_program_id_to_input_id;//v[cell_program_id]=cell_input_id
	unordered_map<cell_input_id_type, cell_program_id_type>cell_input_id_to_program_id;//v[cell_input_id]=cell_program_id
	vector<vector<cell_program_id_type>> netlist;//netlist[0]  for n1,v[net_id]={cell_program_id:cell_program_id conneted net_id}
	vector< vector<net_id_type>>netlist_T;//v[cell_program_id]={net_id:net_id conneted cell_program_id}
};
myGraph file_read_g(const string cell_file, const string net_file)
{
	//genTestcase.rb maybe generate a  cell that no net connects it,so we need to check
	myGraph G;
	ifstream fin;
	fin.open(cell_file);
	string cell_id, area;
	unordered_map<cell_input_id_type, int>cell_input_id_to_area;//v[cell_input_id]=area
	while (fin >> cell_id)//O(#cells)
	{
		fin >> area;
		cell_input_id_to_area[stoi(cell_id.substr(1))] = stoi(area);
	}
	int i = 0;
	for (auto [cell_input_id, area] : cell_input_id_to_area) {
		G.cell_program_id_to_area.push_back(area);
		G.cell_program_id_to_input_id.push_back(cell_input_id);
		G.cell_input_id_to_program_id[cell_input_id] = i;//O(1)
		G.netlist_T.push_back({});
		i++;
	}
	fin.close();
	fin.open(net_file);
	string t;
	while (fin >> t)
	{
		fin >> t >> t;
		vector<cell_program_id_type> cells;
		while (fin >> t)
		{
			if (t == "}")break;
			cell_program_id_type cell_id = G.cell_input_id_to_program_id[stoi(t.substr(1))];
			cells.push_back(cell_id);
			G.netlist_T[cell_id].push_back(G.netlist.size());
		}
		G.netlist.push_back(move(cells));
	}
	fin.close();
	return G;
}
int file_read_partition(vector< partition_name_type>& cell_id_to_partition_name, const string partition_file, const
	myGraph& G) {
	string t;
	ifstream fin;
	fin.open(partition_file);
	fin >> t >> t;
	int cut_size = stoi(t);
	cout << "initial cost:" << cut_size << endl;
	fin >> t >> t;
	for (int num = stoi(t); num > 0; num--) {
		fin >> t;
		int cell_id = G.cell_input_id_to_program_id.at(stoi(t.substr(1)));//O(1)
		cell_id_to_partition_name[cell_id] = A;
	}
	fin >> t >> t;
	for (int num = stoi(t); num > 0; num--) {
		fin >> t;
		int cell_id = G.cell_input_id_to_program_id.at(stoi(t.substr(1)));//O(1)
		cell_id_to_partition_name[cell_id] = B;
	}
	fin.close();
	return  cut_size;
}
pair<int, int> my_BALANCE_CRITERION(const vector<int>& cell_program_id_to_area, const double r, const double epsilon) {//O(#cells)
	int sum = 0, max_area = 0;
	for (const int area : cell_program_id_to_area) {
		sum += area;
		if (area > max_area) max_area = area;
	}
	int lb = min(floor(r * sum) - max_area, floor((r - epsilon) * sum));//min(0.48*sum_area,0.5*sum_area-max_area)
	return { lb, sum - lb };
}
array<int, 2>cell_id_to_partition_name_generate_partition_name_to_partition_area(
	const vector<partition_name_type>& cell_id_to_partition_name, const vector<int>& cell_program_id_to_area) {
	std::array<int, 2>partition_name_to_partition_area = { 0,0 };
	for (int cell_id = 0; cell_id < cell_id_to_partition_name.size(); cell_id++) {
		partition_name_to_partition_area[cell_id_to_partition_name[cell_id]] += cell_program_id_to_area[cell_id];
	}
	return partition_name_to_partition_area;
}
void cell_id_to_partition_name_generate_net_id_to_partition_name_to_cell_num(vector<array< int, 2>>& net_id_to_partition_name_to_cell_num,
	const vector< partition_name_type>& cell_id_to_partition_name, const vector<vector<cell_program_id_type>>& netlist) {
	for (int i = 0; i < netlist.size(); i++) {
		net_id_to_partition_name_to_cell_num[i] = { 0,0 };
		for (int cell_id : netlist[i])
			net_id_to_partition_name_to_cell_num[i][cell_id_to_partition_name[cell_id]]++;
	}
}
vector<int>get_random_i1_and_no_repeat_vector_from_i2(const int i1, vector<int> i2, IRandomIntGenerator& randGen) {
	vector<int>ret(i1);
	if (i1 > i2.size()) {
		cerr << "Error: i1 is larger than the size of i2." << endl;
		exit(EXIT_FAILURE);
	}
	if (i1 > i2.size() / 2) {
		for (int i = 0, i2_size_1 = i2.size() - 1; i < i1; i++) {
			int random_index = randGen.Next(i, i2_size_1);
			ret[i] = i2[random_index];
			i2[random_index] = i2[i];
		}
	}
	else {
		unordered_set<int> selected_int;
		for (int i = 0; i < i1; i++) {
			int random_index;
			do {
				random_index = randGen.Next(0, i2.size() - 1);
			} while (selected_int.count(random_index) > 0);//如果已經選過了，重新抽取//碰撞機率<0.5
			selected_int.insert(random_index);
			ret[i] = i2[random_index];
		}
	}
	return ret;
}
pair<int,int> compute_cut_size_and_A_area_if_move(array< int, 2> partition_name_to_partition_area,
	vector<array< int, 2>>net_id_to_partition_name_to_cell_num, const vector<partition_name_type>& cell_id_to_partition_name, const vector<cell_program_id_type>& V2, const int
	move_int, const myGraph& G, const int lb, const int ub) {
	int pow_2 = 1;
	for (const int cell_id : V2) {
		if (move_int & pow_2) {
			partition_name_to_partition_area[cell_id_to_partition_name[cell_id]] -= G.cell_program_id_to_area[cell_id];
			partition_name_to_partition_area[1 - cell_id_to_partition_name[cell_id]] += G.cell_program_id_to_area[cell_id];

		}
		pow_2 <<= 1;
	}
	if (partition_name_to_partition_area[A] > ub || partition_name_to_partition_area[A] < lb)
		return { INT_MAX,INT_MAX };//不符合平衡條件，返回MAX
	pow_2 = 1;
	for (const int cell_id : V2) {
		if (move_int & pow_2) {
			const partition_name_type& from = cell_id_to_partition_name[cell_id], to = 1 - from;
			for (const int net : G.netlist_T[cell_id]) {//O(degree(cell_id))
				auto& partition_name_to_cell_num_of_net = net_id_to_partition_name_to_cell_num[net];
				partition_name_to_cell_num_of_net[from]--;
				partition_name_to_cell_num_of_net[to]++;
			}
		}
		pow_2 <<= 1;
	}
	int cut_size = 0;
	for (auto& partition_name_to_cell_num_of_net : net_id_to_partition_name_to_cell_num)//O(#nets)
		if (partition_name_to_cell_num_of_net[A] > 0 && partition_name_to_cell_num_of_net[B] > 0)
			cut_size++;
	return { cut_size, labs(partition_name_to_partition_area[A] - partition_name_to_partition_area[B]) };
}
void MOVE(vector<partition_name_type>& cell_id_to_partition_name, const
	vector<cell_program_id_type>& V2, const int move_int, const myGraph& G) {
	for (int m = 0, pow_2_m = 1; m < V2.size(); m++, pow_2_m <<= 1) {
		if (move_int & pow_2_m) {
			cell_program_id_type cell_id = V2[m];
			const partition_name_type from = cell_id_to_partition_name[cell_id], to = A + B - from;
			cell_id_to_partition_name[cell_id] = to;
		}
	}
	return;
}
int main()
{
	unique_ptr<IRandomIntGenerator> randGen_p = std::make_unique<RandModGenerator>(123);
	ofstream delta_cut_size_out("delta_cut_size.txt");
	for (int test_case = 1; test_case <= TEST_CASE_NUM; test_case++) {
		cout << "Test case " << test_case << "\n";
		string file_name_base = /*string("../../colab_testcase 10000 100 10000 2/test") + to_string(test_case);
			//*/"../../ISPD_benchmark/ibm01";//"../../../cpp/src/src/ISPD_benchmark/ibm01";
		string cell_file = file_name_base + ".cells";
		string net_file = file_name_base + ".nets";
		string partition_file = file_name_base + ".partitions";
		const myGraph G = file_read_g(cell_file, net_file);
		const vector<int>range = [](int n)->vector<int>
			{vector<int> V1(n); iota(V1.begin(), V1.end(), 0); return V1; }(G.cell_program_id_to_area.size());
		vector< partition_name_type> cell_id_to_partition_name(G.cell_program_id_to_area.size());
		const int initial_cut_size = file_read_partition(cell_id_to_partition_name, partition_file, G);
		int total_area = 0;
		for (int area : G.cell_program_id_to_area)
			total_area += area;
		const int ub = floor(0.55 * total_area), lb = total_area - ub;//definition of github.com/EricLu1218/Physical_Design_Automation/blob/main/Two-way_Min-cut_Partitioning/CS613500_HW2_spec.pdf
		//const auto [lb, ub] = my_BALANCE_CRITERION(G.cell_program_id_to_area, 0.5, 0.02);//c++17 auto[]
		int min_cut_size_for_testcase = initial_cut_size;
		vector<array< int, 2>>net_id_to_partition_name_to_cell_num(G.netlist.size(), { 0,0 });
		const int round_num = 10;
		double total_time = 0;
		const int enumerate_cell_num = min(max_enumerate_cell_num, (int)G.cell_program_id_to_area.size());
		const int move_method_num = 1 << enumerate_cell_num;
		for (unsigned round = 1; round <= round_num; round++) {
			cout << "Round " << round << "\n";
			array< int, 2>partition_name_to_partition_area =
				cell_id_to_partition_name_generate_partition_name_to_partition_area(cell_id_to_partition_name, G.cell_program_id_to_area);//O(#cells)
			cell_id_to_partition_name_generate_net_id_to_partition_name_to_cell_num(net_id_to_partition_name_to_cell_num, cell_id_to_partition_name, G.netlist);//O(maxDegree * #nets)
			vector< cell_program_id_type>enumerate_cells = get_random_i1_and_no_repeat_vector_from_i2(enumerate_cell_num, range, *randGen_p);

			auto start = chrono::high_resolution_clock::now();
			int local_min_cut_cost = INT_MAX, local_min_imbalance_cost = INT_MAX, local_second_min_cut_cost = INT_MAX, 
				global_min_cut_cost = INT_MAX, global_min_imbalance_cost = INT_MAX, global_second_min_cut_cost = INT_MAX;
			int local_min_idx = -1, global_min_idx = -1, local_second_min_idx = -1, global_second_min_idx = -1;
//#pragma omp parallel /**/num_threads(4) firstprivate(local_min_cut_cost,local_second_min_cut_cost,local_min_imbalance_cost,local_min_idx)
			{
//#pragma omp for schedule(dynamic, 10000)
				for (int i = 0; i < move_method_num; i++) {//在lb,ub外的case比在lb,ub內的case工作量短很多
					auto[cut_cost,imbalance_cost] = compute_cut_size_and_A_area_if_move(partition_name_to_partition_area,
						net_id_to_partition_name_to_cell_num, cell_id_to_partition_name, enumerate_cells, i, G, lb, ub);
					
					if (cut_cost < local_min_cut_cost) {
						cout << "local_second_min_cut_cost change:" << local_second_min_cut_cost << "->" << local_min_cut_cost << "\n";
						local_second_min_cut_cost = local_min_cut_cost;
						local_second_min_idx = local_min_idx;
						local_min_cut_cost = cut_cost;
						local_min_imbalance_cost = imbalance_cost;
						local_min_idx =i;
					}
					else if (cut_cost == local_min_cut_cost ){
						if (imbalance_cost < local_min_imbalance_cost) {
							local_min_imbalance_cost = imbalance_cost;
							local_min_idx = i;
						}
					}
					else if (cut_cost < local_second_min_cut_cost) {
						cout << "local_second_min_cut_cost change:" << local_second_min_cut_cost << "->" << cut_cost << "\n";
						local_second_min_cut_cost = cut_cost;
						local_second_min_idx = i;
					}
				}
//#pragma omp critical
				{
					if (local_min_cut_cost < global_min_cut_cost) {
						cout << "global_min_cut_cost change:" << global_min_cut_cost << "->" << local_min_cut_cost << "\n";
						//global_second_min_cut_cost = min(global_min_cut_cost, local_second_min_cut_cost);
						if(local_second_min_cut_cost<global_min_cut_cost){
							global_second_min_cut_cost=local_second_min_cut_cost;
							global_second_min_idx = local_second_min_idx;
						}
						else {
							global_second_min_cut_cost = global_min_cut_cost;
							global_second_min_idx = global_min_idx;
						}
						global_min_cut_cost = local_min_cut_cost;
						global_min_imbalance_cost = local_min_imbalance_cost;
						global_min_idx = local_min_idx;
					}
					else if (local_min_cut_cost == global_min_cut_cost ){
						if (local_min_imbalance_cost < global_min_imbalance_cost) {
							global_min_imbalance_cost = local_min_imbalance_cost;
							global_min_idx = local_min_idx;
						}
					}
					else if (local_min_cut_cost < global_second_min_cut_cost) {
						global_second_min_cut_cost = local_min_cut_cost;
						global_second_min_idx = local_min_idx;
					}
				}
			}
			min_cut_size_for_testcase = floor(global_min_cut_cost);//建議改用pair<int,int>{cut_cost,unbalance_score}，因為double可能會有誤差
			MOVE(cell_id_to_partition_name, enumerate_cells, global_min_idx, G);
			auto end = chrono::high_resolution_clock::now();
			auto duration = chrono::duration_cast<chrono::milliseconds>(end - start);
			cout << "round:" << round << ".Time taken: " << duration.count() << " ms" << "\n";
			total_time += duration.count();
			vector<int> area(2);
			for (int i = 0; i < cell_id_to_partition_name.size(); i++)
				area[cell_id_to_partition_name[i]] += G.cell_program_id_to_area[i];
			cout << "A area:" << area[A] << "\n" << "B area:" << area[B] << "\n";
			cout << "test_min_cut_size:" << min_cut_size_for_testcase << endl 
				<< "test_min_idx:" << global_min_idx << "\n"
				<< "test_second_min_cut_size:" << global_second_min_cut_cost << "\n"
				<< "test_second_min_idx:" << global_second_min_idx << "\n";
		}
		cout << "Average time per round: " << total_time / round_num << " ms\nmin cut size:" << min_cut_size_for_testcase << "\n";
		delta_cut_size_out << initial_cut_size << "->" << min_cut_size_for_testcase << "\n";
	}
	delta_cut_size_out.close();
	return 0;
}