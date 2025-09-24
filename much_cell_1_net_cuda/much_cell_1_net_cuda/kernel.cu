using namespace std;
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
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
#include <execution>
constexpr int max_enumerate_cell_num = 20;
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
	vector<int>cell_program_id_to_area;//a[cell_program_id]=area
	vector<cell_input_id_type>cell_program_id_to_input_id;//a[cell_program_id]=cell_input_id
	unordered_map<cell_input_id_type, cell_program_id_type>cell_input_id_to_program_id;//v[cell_input_id]=cell_program_id
	vector<vector<cell_program_id_type>> netlist;//netlist[0]  for n1,v[net_id]={cell_program_id:cell_program_id conneted net_id}
	vector< vector<net_id_type>>netlist_T;//v[cell_program_id]={net_id:net_id conneted cell_program_id}
	int total_area, total_edge_weight;
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
	for (auto &a : cell_input_id_to_area) {
		auto cell_input_id = a.first, area = a.second;
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
	G.total_area = accumulate(G.cell_program_id_to_area.begin(), G.cell_program_id_to_area.end(), 0);
	G.total_edge_weight = G.netlist.size();//每條邊權重1
	return G;
}
int file_read_partition (vector< partition_name_type> * cell_id_to_partition_name,const string& partition_file, const
myGraph& G) {
	cell_id_to_partition_name->resize(G.cell_program_id_to_area.size());
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
		(*cell_id_to_partition_name)[cell_id] = A;
	}
	fin >> t >> t;
	for (int num = stoi(t); num > 0; num--) {
		fin >> t;
		int cell_id = G.cell_input_id_to_program_id.at(stoi(t.substr(1)));//O(1)
		(*cell_id_to_partition_name)[cell_id] = B;
	}
	fin.close();
	return  cut_size ;
}
pair<int, int> my_BALANCE_CRITERION(const myGraph& G, const /*double r*/unsigned P , const /*double epsilon*/unsigned x) {//O(#cells)
	//definition of https://chriswalshaw.co.uk/partition/
	int  S_opt = (G.total_area + P - 1) / P;//ceil(1.0*graph.total_area/P)用double不准
	int max_S_p = S_opt * (100 + x) / 100;
	return { G.total_area - max_S_p, max_S_p };
}
array<int, 2>cell_id_to_partition_name_generate_partition_name_to_partition_area(
	const vector<partition_name_type>& cell_id_to_partition_name, const vector<int>& cell_program_id_to_area) {
	std::array<int, 2>partition_name_to_partition_area = { 0,0 };
	for (int cell_id = 0; cell_id < cell_id_to_partition_name.size(); cell_id++) {
		partition_name_to_partition_area[cell_id_to_partition_name[cell_id]] += cell_program_id_to_area[cell_id];
	}
	return partition_name_to_partition_area;
}
vector<array< int, 2>> cell_id_to_partition_name_generate_net_id_to_partition_name_to_cell_num(
	const vector< partition_name_type>& cell_id_to_partition_name, const myGraph& G) {
	vector<array< int, 2>> net_id_to_partition_name_to_cell_num(G.netlist.size(), { 0,0 });
	for (int i = 0; i < G.netlist.size(); i++) {
		net_id_to_partition_name_to_cell_num[i] = { 0,0 };
		for (int cell_id :G.netlist[i])
			net_id_to_partition_name_to_cell_num[i][cell_id_to_partition_name[cell_id]]++;
	}
	return move(net_id_to_partition_name_to_cell_num);
}
vector<int>get_random_i1_and_no_repeat_vector_from_0_to_i2(//const myGraph& G,
	const int num, int max_plus_1, IRandomIntGenerator& randGen) {
	const static vector<int>range = [](int n)->vector<int> {vector<int> V1(n); iota(V1.begin(), V1.end(), 0); return V1; }(max_plus_1);

	vector<int>ret(num);
	if (num > max_plus_1) {
		cerr << "Error: i1 is larger than the size of i2." << endl;
		exit(EXIT_FAILURE);
	}
	if (num > max_plus_1 / 2) {
		auto range_copy = range;
		for (int i = 0, i2_size_1 = max_plus_1 - 1; i < num; i++) {
			int random_index = randGen.Next(i, i2_size_1);
			ret[i] = range_copy[random_index];
			range_copy[random_index] = range_copy[i];
		}
	}
	else {
		unordered_set<int> selected_int;
		for (int i = 2; i < num; i++) {
			int random_index;
			do {
				random_index = randGen.Next(0, max_plus_1 - 1);
			} while (selected_int.count(random_index) > 0);//如果已經選過了，重新抽取//碰撞機率<0.5
			selected_int.insert(random_index);
			ret[i] = range[random_index];
		}
	}
	return ret;
}
pair<int, int> compute_cut_size_and_A_area_if_move(array< int, 2> partition_name_to_partition_area,
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
struct cost_type {
	int cut_cost;
	int balance_cost;
	unsigned move_idx;
	__host__ __device__ bool operator<(const cost_type& right) {
		if (cut_cost < right.cut_cost)return true;
		if (cut_cost == right.cut_cost && balance_cost < right.balance_cost)return true;
		return false;
	}
};
int main()
{
	unique_ptr<IRandomIntGenerator> randGen_p = std::make_unique<RandModGenerator>(123);
	ofstream delta_cut_size_out("delta_cut_size.txt");
	//for (int test_case = 1; test_case <= TEST_CASE_NUM; test_case++) {
	//cout << "Test case " << test_case << "\n";
	string file_name_base = /*string("../../colab_testcase 10000 100 10000 2/test") + to_string(test_case);
		//*/"../../ISPD_benchmark/ibm01";
	string cell_file = file_name_base + ".cells";
	string net_file = file_name_base + ".nets";
	string partition_file = file_name_base + ".partitions";
	const myGraph G = file_read_g(cell_file, net_file);
	vector< partition_name_type> cell_id_to_partition_name; 
	int initial_cut_size= file_read_partition(&cell_id_to_partition_name,partition_file, G);
	int lb, ub;
	tie(lb, ub) = my_BALANCE_CRITERION(G, 2, 2);////definition of https://chriswalshaw.co.uk/partition/
	vector<array< int, 2>>net_id_to_partition_name_to_cell_num=
		cell_id_to_partition_name_generate_net_id_to_partition_name_to_cell_num(  cell_id_to_partition_name, G);
	cout << "initial cut size:" << initial_cut_size << "\n";
	int min_cut_size_for_testcase = initial_cut_size;
	const int round_num = 100;
	double total_time = 0;
	const int enumerate_cell_num = min(max_enumerate_cell_num, (int)G.cell_program_id_to_area.size());
	const int move_method_num = 1 << enumerate_cell_num;
	for (unsigned round = 1; round <= round_num; round++) {
		cout << "Round " << round << "\n";
		array< int, 2>partition_name_to_partition_area =
			cell_id_to_partition_name_generate_partition_name_to_partition_area(cell_id_to_partition_name, G.cell_program_id_to_area);//O(#cells)
		net_id_to_partition_name_to_cell_num=cell_id_to_partition_name_generate_net_id_to_partition_name_to_cell_num(
			cell_id_to_partition_name, G);//O(maxDegree * #nets)
		vector< cell_program_id_type>enumerate_cells = get_random_i1_and_no_repeat_vector_from_0_to_i2(
			enumerate_cell_num, G.cell_program_id_to_area.size(), *randGen_p);
		array<int,2> affected_part_area=[&]() {
			array<int, 2> affected_part_area = { 0,0 };
			for (const int cell_id : enumerate_cells) {
				affected_part_area[cell_id_to_partition_name[cell_id]] += G.cell_program_id_to_area[cell_id];
			}
			return affected_part_area;
			}();
		vector<net_id_type> affected_nets=[&]() {
			unordered_set<net_id_type> affected_nets_set;
			for (const int cell_id : enumerate_cells) {
				for (const int net_id : G.netlist_T[cell_id]) {
					affected_nets_set.insert(net_id);
				}
			}
			return vector<net_id_type>(affected_nets_set.begin(), affected_nets_set.end());
			}();
		int affected_cut_size = [&]() {
			int cut_size = 0;
			for (const int net_id : affected_nets) {
				if (net_id_to_partition_name_to_cell_num[net_id][A] > 0 && net_id_to_partition_name_to_cell_num[net_id][B] > 0)
					cut_size++;
			}
			return cut_size;
			}();
		auto start = chrono::high_resolution_clock::now();
		cost_type global_min_cost = {G.total_edge_weight,G.total_area,INT_MAX}, global_second_min_cost= global_min_cost;
		//#pragma omp parallel /**/num_threads(4) firstprivate(local_min_cut_cost,local_second_min_cut_cost,local_min_imbalance_cost,local_min_idx)
		{
			//#pragma omp for schedule(dynamic, 10000)
			for (int i = 0; i < move_method_num; i++) {//在lb,ub外的case比在lb,ub內的case工作量短很多
				auto [cut_cost, imbalance_cost] = compute_cut_size_and_A_area_if_move(partition_name_to_partition_area,
					net_id_to_partition_name_to_cell_num, cell_id_to_partition_name, enumerate_cells, i, G, lb, ub);

				if (cut_cost < local_min_cut_cost) {
					cout << "local_second_min_cut_cost change:" << local_second_min_cut_cost << "->" << local_min_cut_cost << "\n";
					local_second_min_cut_cost = local_min_cut_cost;
					local_second_min_idx = local_min_idx;
					local_min_cut_cost = cut_cost;
					local_min_imbalance_cost = imbalance_cost;
					local_min_idx = i;
				}
				else if (cut_cost == local_min_cut_cost) {
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
					if (local_second_min_cut_cost < global_min_cut_cost) {
						global_second_min_cut_cost = local_second_min_cut_cost;
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
				else if (local_min_cut_cost == global_min_cut_cost) {
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
	//}
	delta_cut_size_out.close();
	return 0;
}
/*
cudaError_t addWithCuda(int *c, const int *a, const int *b, unsigned int size);

__global__ void addKernel(int *c, const int *a, const int *b)
{
    int i = threadIdx.x;
    c[i] = a[i] + b[i];
}

int main()
{
    const int arraySize = 5;
    const int a[arraySize] = { 1, 2, 3, 4, 5 };
    const int b[arraySize] = { 10, 20, 30, 40, 50 };
    int c[arraySize] = { 0 };

    // Add vectors in parallel.
    cudaError_t cudaStatus = addWithCuda(c, a, b, arraySize);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "addWithCuda failed!");
        return 1;
    }

    printf("{1,2,3,4,5} + {10,20,30,40,50} = {%d,%d,%d,%d,%d}\n",
        c[0], c[1], c[2], c[3], c[4]);

    // cudaDeviceReset must be called before exiting in order for profiling and
    // tracing tools such as Nsight and Visual Profiler to show complete traces.
    cudaStatus = cudaDeviceReset();
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaDeviceReset failed!");
        return 1;
    }

    return 0;
}

// Helper function for using CUDA to add vectors in parallel.
cudaError_t addWithCuda(int *c, const int *a, const int *b, unsigned int size)
{
    int *dev_a = 0;
    int *dev_b = 0;
    int *dev_c = 0;
    cudaError_t cudaStatus;

    // Choose which GPU to run on, change this on a multi-GPU system.
    cudaStatus = cudaSetDevice(0);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaSetDevice failed!  Do you have a CUDA-capable GPU installed?");
        goto Error;
    }

    // Allocate GPU buffers for three vectors (two input, one output)    .
    cudaStatus = cudaMalloc((void**)&dev_c, size * sizeof(int));
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMalloc failed!");
        goto Error;
    }

    cudaStatus = cudaMalloc((void**)&dev_a, size * sizeof(int));
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMalloc failed!");
        goto Error;
    }

    cudaStatus = cudaMalloc((void**)&dev_b, size * sizeof(int));
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMalloc failed!");
        goto Error;
    }

    // Copy input vectors from host memory to GPU buffers.
    cudaStatus = cudaMemcpy(dev_a, a, size * sizeof(int), cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMemcpy failed!");
        goto Error;
    }

    cudaStatus = cudaMemcpy(dev_b, b, size * sizeof(int), cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMemcpy failed!");
        goto Error;
    }

    // Launch a kernel on the GPU with one thread for each element.
    addKernel<<<1, size>>>(dev_c, dev_a, dev_b);

    // Check for any errors launching the kernel
    cudaStatus = cudaGetLastError();
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "addKernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
        goto Error;
    }
    
    // cudaDeviceSynchronize waits for the kernel to finish, and returns
    // any errors encountered during the launch.
    cudaStatus = cudaDeviceSynchronize();
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching addKernel!\n", cudaStatus);
        goto Error;
    }

    // Copy output vector from GPU buffer to host memory.
    cudaStatus = cudaMemcpy(c, dev_c, size * sizeof(int), cudaMemcpyDeviceToHost);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMemcpy failed!");
        goto Error;
    }

Error:
    cudaFree(dev_c);
    cudaFree(dev_a);
    cudaFree(dev_b);
    
    return cudaStatus;
}
*/