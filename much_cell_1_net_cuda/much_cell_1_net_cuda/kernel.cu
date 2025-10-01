using namespace std;
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <iostream>
#include<fstream>
#include<numeric>
#include <cfloat>
#include<tuple>
#include <vector>
#include <string>
#include<climits>
#include<unordered_set>
#include<unordered_map>
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
	vector<int>cellId_to_area;//a[cell_program_id]=area
	vector<int>cellId_to_cellInputId;//a[cell_program_id]=cell_input_id
	unordered_map<int, int>cellInputId_to_cellId;//a[cell_input_id]=cell_program_id
	vector<vector<int>> netId_to_cellId_s;//[0]  for n1,a[net_id]={cell_program_id:cell_program_id conneted net_id}
	vector< vector<int>>cellId_to_netId_s;//a[cell_program_id]={net_id:net_id conneted cell_program_id}
	int total_area, total_edge_weight;
};
myGraph file_read_g(const string cell_file, const string net_file)
{
	//genTestcase.rb maybe generate a  cell that no net connects it,so we need to check
	myGraph G;
	ifstream fin;
	fin.open(cell_file);
	string cell_id, area;
	unordered_map<int, int>cellId_to_area;//a[cell_input_id]=area
	while (fin >> cell_id)//O(#cells)
	{
		fin >> area;
		cellId_to_area[stoi(cell_id.substr(1))] = stoi(area);
	}
	int i = 0;
	for (auto& a : cellId_to_area) {
		auto cell_input_id = a.first, area = a.second;
		G.cellId_to_area.push_back(area);
		G.cellId_to_cellInputId.push_back(cell_input_id);
		G.cellInputId_to_cellId[cell_input_id] = i;//O(1)
		G.cellId_to_netId_s.push_back({});
		i++;
	}
	fin.close();
	fin.open(net_file);
	string t;
	while (fin >> t)
	{
		fin >> t >> t;
		vector<int> cellId_s;
		while (fin >> t)
		{
			if (t == "}")break;
			int cell_id = G.cellInputId_to_cellId[stoi(t.substr(1))];
			cellId_s.push_back(cell_id);
			G.cellId_to_netId_s[cell_id].push_back(G.netId_to_cellId_s.size());
		}
		G.netId_to_cellId_s.push_back(move(cellId_s));
	}
	fin.close();
	G.total_area = accumulate(G.cellId_to_area.begin(), G.cellId_to_area.end(), 0);
	G.total_edge_weight = G.netId_to_cellId_s.size();//每條邊權重1
	return G;
}
int file_read_partition(vector< int>* cellId_to_partId, const string& partition_file, const
	myGraph& G) {
	cellId_to_partId->resize(G.cellId_to_area.size());
	string t;
	ifstream fin;
	fin.open(partition_file);
	fin >> t >> t;
	int cut_size = stoi(t);
	cout << "initial cost:" << cut_size << endl;
	fin >> t >> t;
	for (int num = stoi(t); num > 0; num--) {
		fin >> t;
		int cell_id = G.cellInputId_to_cellId.at(stoi(t.substr(1)));//O(1)
		(*cellId_to_partId)[cell_id] = A;
	}
	fin >> t >> t;
	for (int num = stoi(t); num > 0; num--) {
		fin >> t;
		int cell_id = G.cellInputId_to_cellId.at(stoi(t.substr(1)));//O(1)
		(*cellId_to_partId)[cell_id] = B;
	}
	fin.close();
	return  cut_size;
}
pair<int, int> my_BALANCE_CRITERION(const myGraph& G, const unsigned P, const unsigned x) {//O(#cells)
	//definition of https://chriswalshaw.co.uk/partition/
	int  S_opt = (G.total_area + P - 1) / P;//ceil(1.0*graph.total_area/P)用double不准
	int max_S_p = S_opt * (100 + x) / 100;
	return { G.total_area - max_S_p, max_S_p };
}
array<int, 2>cell_id_to_partition_name_generate_partition_name_to_partition_area(
	const vector<int>& cellId_to_partId, const vector<int>& cellId_to_area) {
	std::array<int, 2>partId_to_partArea = { 0,0 };
	for (int cell_id = 0; cell_id < cellId_to_partId.size(); cell_id++) {
		partId_to_partArea[cellId_to_partId[cell_id]] += cellId_to_area[cell_id];
	}
	return partId_to_partArea;
}
vector<array< int, 2>> cell_id_to_partition_name_generate_net_id_to_partition_name_to_cell_num(
	const vector< int>& cellId_to_partId, const myGraph& G) {
	vector<array< int, 2>> netId_to_partId_to_cellNum(G.netId_to_cellId_s.size(), { 0,0 });
	for (int i = 0; i < G.netId_to_cellId_s.size(); i++) {
		netId_to_partId_to_cellNum[i] = { 0,0 };
		for (int cell_id : G.netId_to_cellId_s[i])
			netId_to_partId_to_cellNum[i][cellId_to_partId[cell_id]]++;
	}
	return move(netId_to_partId_to_cellNum);
}
vector<int>get_random_i1_and_no_repeat_vector_from_0_to_i2_minus_1(//const myGraph& G,
	const int num, int max_plus_1, IRandomIntGenerator& randGen) {
	const static vector<int>range = [](int n)->vector<int> {vector<int> V1(n); iota(V1.begin(), V1.end(), 0); return V1; }(max_plus_1);

	vector<int>ret(num);
	if (num > max_plus_1) {
		cerr << "Error: i1 is larger than the size of i2." << endl;
		exit(EXIT_FAILURE);
	}
	if (num > max_plus_1 / 2) {
		vector<int> range_copy = range;
		for (int i = 0, i2_size_1 = max_plus_1 - 1; i < num; i++) {
			int random_index = randGen.Next(i, i2_size_1);
			ret[i] = range_copy[random_index];
			range_copy[random_index] = range_copy[i];
		}
	}
	else {
		unordered_set<int> selected_int;
		for (int i = 0; i < num; i++) {
			int random_index;
			do {
				random_index = randGen.Next(0, max_plus_1 - 1);
			} while (selected_int.count(random_index) > 0);//如果已經選過了，重新抽取//碰撞機率<0.5
			selected_int.insert(random_index);
			ret[i] = range[random_index];
		}
	}
	return move(ret);
}
void MOVE(vector<int>* cell_id_to_partition_name, const
	vector<int>& enumerate_cellId_s, const int move_int, const myGraph& G) {
	for (int m = 0, pow_2_m = 1; m < enumerate_cellId_s.size(); m++, pow_2_m <<= 1) {
		/*if (move_int & pow_2_m) {
			int cell_id = enumerate_cellId_s[m];
			const int from = (*cell_id_to_partition_name)[cell_id], to = A + B - from;
			(*cell_id_to_partition_name)[cell_id] = to;
		}*/
		int cell_id = enumerate_cellId_s[m];
		(*cell_id_to_partition_name)[cell_id] =(int)( (move_int & pow_2_m)!=0);
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
template<typename T>
struct myVector {
	T*data;
	int size;
	__host__ __device__ T& operator[](int idx) { return data[idx]; }
	__host__ __device__ T at(int idx) const { return data[idx]; }
};
struct device_graph {
	myVector<int>cellId_to_area;
	myVector<myVector<int>>netId_to_cellId_s;
	myVector<myVector<int>>cellId_to_netId_s;
	int total_area, total_edge_weight;
};
struct device_part {
	myVector<int>cellId_to_partId;
	int partId_to_area[2];
	int ub;
};
struct device_task {
	myVector<int>nodes;
	myVector<int>edges;
	struct hash_node {
		int node_idx, task_idx, next_idx;
	}*task_node_at_i_to_i_data;
	__device__ int nodeIdx_to_taskIdx(int node_idx) const {
		int hash_idx = node_idx % (2 * nodes.size);
		//linked list,data[0~2*nodes_num-1} are dummy heads,data[2*nodes_num~] are list nodes,load factor=1/2
		hash_idx = task_node_at_i_to_i_data[hash_idx].next_idx;
		while (hash_idx != -1) {
			if (task_node_at_i_to_i_data[hash_idx].node_idx == node_idx)return task_node_at_i_to_i_data[hash_idx].task_idx;
			hash_idx = task_node_at_i_to_i_data[hash_idx].next_idx;
		}
		return nodes.size;
	}
};
__device__ cost_type device_evaluate(const device_graph& graph,const device_part& extracted_part,
	const device_task&task, unsigned move_idx,unsigned not_move_idx) {
	//if (move_idx == not_move_idx)printf("task nodes size %d, edges size %d\n", task.nodes.size, task.edges.size);
	unsigned area[2] = { 0, 0 };
	for (unsigned i = 0; i < task.nodes.size; i++) {
		unsigned part = ((move_idx & (1 << i)) != 0);
		area[part] += graph.cellId_to_area.at(task.nodes.at(i));
		/*for (int i = 0, size = graph.cellId_to_netId_s.at(task.nodes.at(i)).size; i < size; i++) {
			auto netId = graph.cellId_to_netId_s.at(task.nodes.at(i)).at(i);
			tmp_netId_to_partId_to_cell_num[netId][part]++;
		}*/
	}
	for (unsigned i = 0; i < 2; i++) 
		if (extracted_part.partId_to_area[i] + area[i] > extracted_part.ub) 
			return  { graph.total_edge_weight ,graph.total_area ,move_idx };
	int balance_cost = labs(extracted_part.partId_to_area[A] + area[A] - extracted_part.partId_to_area[B] - area[B]);
	int cut_cost = 0;
	//if (move_idx == not_move_idx)printf("calculate not_move_idx %d \n", not_move_idx);
	for (int j = 0; j< task.edges.size; j++) {
		unsigned edgeId = task.edges.at(j);
		int temp_partIdx_to_cell_num[2] = { 0,0 };
		int idx;
		for (int i = 0, size = graph.netId_to_cellId_s.at(edgeId).size; i < size; i++) {
			int cellId = graph.netId_to_cellId_s.at(edgeId)[i];
			unsigned part = ((idx = task.nodeIdx_to_taskIdx(cellId)) < task.nodes.size)
				? ((move_idx & (1 << idx)) != 0) :extracted_part.cellId_to_partId.at(cellId);
			temp_partIdx_to_cell_num[part]++;
		}
		if (temp_partIdx_to_cell_num[0] > 0 && temp_partIdx_to_cell_num[1] > 0){
			cut_cost++;//edge weight=1
		}
			//if (move_idx == not_move_idx)printf("%d:%d,%d\n",  edgeId, temp_partIdx_to_cell_num[0], temp_partIdx_to_cell_num[1]);
	}
	//if (move_idx == not_move_idx)printf("\nedges num:%d\n",task.edges.size);
	return  { cut_cost, balance_cost,move_idx };
}
__global__  void evaluate_array_kernel(cost_type* cost, device_part* dev_extracted_part_p, device_graph* device_graph1_p,
	device_task* device_task1_p,  int num_explores,unsigned not_move_idx) {
	unsigned cost_idx = blockIdx.x * blockDim.x + threadIdx.x, move_idx = cost_idx;
	unsigned stride = blockDim.x * gridDim.x;
	cost[cost_idx] = { device_graph1_p->total_edge_weight,device_graph1_p->total_area,0 };
	for (; move_idx < num_explores; move_idx += stride) {
		cost_type now_cost = device_evaluate(*device_graph1_p,*dev_extracted_part_p, *device_task1_p, move_idx, not_move_idx);
		//if (move_idx == not_move_idx)printf("thread %d calculate not_move_idx %d cost %d,%d\n", cost_idx, not_move_idx, now_cost.cut_cost, now_cost.balance_cost);
		if (now_cost < cost[cost_idx])cost[cost_idx] = now_cost;
	}
}
int main()
{
	cudaError_t cudaStatus = cudaSetDevice(0);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaSetDevice failed!  Do you have a CUDA-capable GPU installed?");
		exit(1);
	}
	unique_ptr<IRandomIntGenerator> randGen_p = std::make_unique<RandModGenerator>(123);
	ofstream cost_fout("cost.txt");
	string file_name_base =/* string("../../colab_testcase 10000 100 10000 10/test1");
		//*/"../../ISPD_benchmark/ibm01";
	string cell_file = file_name_base + ".cells";
	string net_file = file_name_base + ".nets";
	string partition_file = file_name_base + ".partitions";
	const myGraph G = file_read_g(cell_file, net_file);
	vector< int> cellId_to_partId;
	/*file_read_partition(&cellId_to_partId, partition_file, G);
	/* //*/
	[&cellId_to_partId, &G]() {
		ifstream part_fin("../../ISPD_benchmark/ibm01.part.2");
		int part_name;
		cellId_to_partId.resize(G.cellId_to_area.size());
		for (int i = 1; i <= G.cellId_to_area.size(); i++) {
			part_fin >> part_name;
			cellId_to_partId[G.cellInputId_to_cellId.at(i)] = part_name;
		}
		part_fin.close();
		}();/**/
	int lb, ub;
	tie(lb, ub) = my_BALANCE_CRITERION(G, 2, 2*2);//definition of https://chriswalshaw.co.uk/partition/
	printf("lb=%d,ub=%d\n", lb, ub);
	vector<array< int, 2>>netId_to_partId_to_cellNum;
	array< int, 2>partId_to_partArea;
	const int round_num = 100;
	double total_time = 0;
	const int enumerate_cell_num = min(max_enumerate_cell_num, (int)G.cellId_to_area.size());
	const int move_method_num = 1 << enumerate_cell_num;
	for (unsigned round = 1; round <= round_num; round++) {
		cout << "Round " << round << "\n";
		partId_to_partArea =
			cell_id_to_partition_name_generate_partition_name_to_partition_area(cellId_to_partId, G.cellId_to_area);		
		cout << "A area = " << partId_to_partArea[A] << " , B area = " << partId_to_partArea[B] << " ,balance_cost="
			<<labs(partId_to_partArea[A]- partId_to_partArea[B])<<"\n";
		netId_to_partId_to_cellNum = cell_id_to_partition_name_generate_net_id_to_partition_name_to_cell_num(
			cellId_to_partId, G);
		int cut_size=[&]() {
			int cut_size = 0;
			for (auto& partition_name_to_cell_num_of_net : netId_to_partId_to_cellNum)//O(#nets)
				if (partition_name_to_cell_num_of_net[A] > 0 && partition_name_to_cell_num_of_net[B] > 0)
					cut_size++;
			return cut_size;
			}();
		cout << "cut size:" << cut_size << "\n";
		vector< int>enumerate_cellId_s = get_random_i1_and_no_repeat_vector_from_0_to_i2_minus_1(
			enumerate_cell_num, G.cellId_to_area.size(), *randGen_p);
		unsigned not_move_idx=[&]() {
			unsigned not_move_idx = 0;
			for (int i = 0; i < enumerate_cellId_s.size(); i++) {
				if (cellId_to_partId[enumerate_cellId_s[i]] == B)
					not_move_idx |= (1 << i);
			}
			return not_move_idx;
			}();
		array<int, 2>extracted_partId_to_area = [&]() {
			array<int, 2>extracted_partId_to_area = { partId_to_partArea[A],partId_to_partArea[B] };
			for (const int cell_id : enumerate_cellId_s) 
					extracted_partId_to_area[cellId_to_partId[cell_id]] -= G.cellId_to_area[cell_id];
			return extracted_partId_to_area;
			}();
		vector<int> affected_netId_s = [&]() {
			unordered_set<int> affected_netId_s_set;
			for (const int cell_id : enumerate_cellId_s) {
				for (const int net_id : G.cellId_to_netId_s[cell_id]) {
					affected_netId_s_set.insert(net_id);
				}
			}
			return vector<int>(affected_netId_s_set.begin(), affected_netId_s_set.end());
			}();
		int affected_cut_size = [&]() {
			int affected_cut_size = 0;
			//cout << "affect_cut_net:";
			for (const int net_id : affected_netId_s) {
				if (netId_to_partId_to_cellNum[net_id][A] > 0 && netId_to_partId_to_cellNum[net_id][B] > 0)
					affected_cut_size++;

				//printf("%d:%d,%d\n", net_id, netId_to_partId_to_cellNum[net_id][A], netId_to_partId_to_cellNum[net_id][B]);
			}
			//cout << "\naffect_net_num:"<< affected_netId_s.size()<<endl;
			//printf("affect cut size %d\n", affected_cut_size);
			return affected_cut_size;
			}();
		auto start = chrono::high_resolution_clock::now();
		vector<cost_type> cost=[&]() {
			int threadsPerBlock = 512,blocksPerGrid = 16;
			cost_type* dev_cost;
			cudaMalloc(&dev_cost, blocksPerGrid * threadsPerBlock * sizeof(cost_type));

			device_task* device_task_p,device_task_cpu;
			cudaMalloc(&device_task_p, sizeof(device_task));

			device_task_cpu.nodes.size = enumerate_cellId_s.size();
			cudaMalloc(&device_task_cpu.nodes.data, enumerate_cellId_s.size() * sizeof(int));
			cudaMemcpy(device_task_cpu.nodes.data, enumerate_cellId_s.data(), enumerate_cellId_s.size() * sizeof(int), cudaMemcpyHostToDevice);
			
			vector<device_task::hash_node> hash_data(3 * enumerate_cellId_s.size());
			for (int i = 0; i < 2 * enumerate_cellId_s.size(); i++)hash_data[i] = device_task::hash_node{ -1,-1,-1 };
			for (int i = 0, free_idx = 2 * enumerate_cellId_s.size(); i < enumerate_cellId_s.size(); i++) {
				int hash_idx = enumerate_cellId_s[i] % (2 * enumerate_cellId_s.size());
				//list head is hash_data[0~2*task.nodes.size()-1], list node is hash_data[2*task.nodes.size()~],load factor=1/2
				hash_data[free_idx] = device_task::hash_node{ enumerate_cellId_s[i], i,hash_data[hash_idx].next_idx };
				hash_data[hash_idx].next_idx = free_idx;
				free_idx++;
			}
			cudaMalloc(&device_task_cpu.task_node_at_i_to_i_data, hash_data.size() * sizeof(device_task::hash_node));
			cudaMemcpy(device_task_cpu.task_node_at_i_to_i_data, hash_data.data(), hash_data.size() * sizeof(device_task::hash_node), cudaMemcpyHostToDevice);
			
			cudaMalloc(&device_task_cpu.edges.data, affected_netId_s.size() * sizeof(int));
			cudaMemcpy(device_task_cpu.edges.data, affected_netId_s.data(), affected_netId_s.size() * sizeof(int), cudaMemcpyHostToDevice);
			device_task_cpu.edges.size = affected_netId_s.size();

			cudaMemcpy(device_task_p, &device_task_cpu, sizeof(device_task), cudaMemcpyHostToDevice);

			device_graph* device_graph_p,device_graph_cpu;
			cudaMalloc(&device_graph_p, sizeof(device_graph));

			device_graph_cpu.total_area = G.total_area;
			device_graph_cpu.total_edge_weight = G.total_edge_weight;

			cudaMalloc(&device_graph_cpu.cellId_to_area.data, G.cellId_to_area.size() * sizeof(int));
			cudaMemcpy(device_graph_cpu.cellId_to_area.data, G.cellId_to_area.data(), G.cellId_to_area.size() * sizeof(int), cudaMemcpyHostToDevice);
			device_graph_cpu.cellId_to_area.size = G.cellId_to_area.size();
			
			device_graph_cpu.cellId_to_netId_s.size = G.cellId_to_netId_s.size();
			cudaMalloc(&device_graph_cpu.cellId_to_netId_s.data, G.cellId_to_netId_s.size() * sizeof(myVector<int>));
			vector<myVector<int>> tmp_cellId_to_netId_s(G.cellId_to_netId_s.size());
			for(int i=0;i<G.cellId_to_netId_s.size();i++){
				auto& tmp = tmp_cellId_to_netId_s[i];
				tmp.size=G.cellId_to_netId_s[i].size();
				cudaMalloc(&tmp.data,tmp.size*sizeof(int));
				cudaMemcpy(tmp.data,G.cellId_to_netId_s[i].data(),tmp.size*sizeof(int),cudaMemcpyHostToDevice);
			}
			cudaMemcpy(device_graph_cpu.cellId_to_netId_s.data, tmp_cellId_to_netId_s.data(), G.cellId_to_netId_s.size() * sizeof(myVector<int>), cudaMemcpyHostToDevice);

			device_graph_cpu.netId_to_cellId_s.size = G.netId_to_cellId_s.size();
			cudaMalloc(&device_graph_cpu.netId_to_cellId_s.data, G.netId_to_cellId_s.size() * sizeof(myVector<int>));
			vector<myVector<int>> tmp_netId_to_cellId_s(G.netId_to_cellId_s.size());
			for(int i=0;i<G.netId_to_cellId_s.size();i++){
				auto& tmp = tmp_netId_to_cellId_s[i];
				tmp.size=G.netId_to_cellId_s[i].size();
				cudaMalloc(&tmp.data,tmp.size*sizeof(int));
				cudaMemcpy(tmp.data,G.netId_to_cellId_s[i].data(),tmp.size*sizeof(int),cudaMemcpyHostToDevice);
				//cudaMemcpy(&device_graph_cpu.netId_to_cellId_s.data[i],&tmp,sizeof(myVector<int>),cudaMemcpyHostToDevice);
			}
			cudaMemcpy(device_graph_cpu.netId_to_cellId_s.data, tmp_netId_to_cellId_s.data(), G.netId_to_cellId_s.size() * sizeof(myVector<int>), cudaMemcpyHostToDevice);

			cudaMemcpy(device_graph_p, &device_graph_cpu, sizeof(device_graph), cudaMemcpyHostToDevice);

			device_part* dev_extracted_part_p,dev_extracted_part_cpu;
			cudaMalloc(&dev_extracted_part_p, sizeof(device_part));

			dev_extracted_part_cpu.ub = ub;
			dev_extracted_part_cpu.partId_to_area[0] = extracted_partId_to_area[0],
				dev_extracted_part_cpu.partId_to_area[1] = extracted_partId_to_area[1];

			dev_extracted_part_cpu.cellId_to_partId.size = cellId_to_partId.size();
			cudaMalloc(&dev_extracted_part_cpu.cellId_to_partId.data, cellId_to_partId.size() * sizeof(int));
			cudaMemcpy(dev_extracted_part_cpu.cellId_to_partId.data, cellId_to_partId.data(), cellId_to_partId.size() * sizeof(int), cudaMemcpyHostToDevice);

			cudaMemcpy(dev_extracted_part_p, &dev_extracted_part_cpu, sizeof(device_part), cudaMemcpyHostToDevice);

			evaluate_array_kernel << <blocksPerGrid, threadsPerBlock >> > (
				dev_cost, dev_extracted_part_p,
				device_graph_p,
				device_task_p,
				move_method_num,not_move_idx
				);

			vector<cost_type> cost;
			cost.resize(blocksPerGrid* threadsPerBlock);
			cudaDeviceSynchronize();
			cudaMemcpy(cost.data(), dev_cost, blocksPerGrid* threadsPerBlock * sizeof(cost_type), cudaMemcpyDeviceToHost);
			/*cout << "not_move_idx:" << not_move_idx << " , mod thread num:" << not_move_idx % (blocksPerGrid *
				threadsPerBlock) << ", cost:" << cost[not_move_idx % (blocksPerGrid *
				threadsPerBlock)].cut_cost + cut_size - affected_cut_size
				<< "," << cost[not_move_idx % (blocksPerGrid *
					threadsPerBlock)].balance_cost << "\n";*/
			cudaFree(dev_cost);

			cudaFree(dev_extracted_part_cpu.cellId_to_partId.data);
			cudaFree(dev_extracted_part_p);

			cudaFree(device_task_cpu.task_node_at_i_to_i_data);
			cudaFree(device_task_cpu.edges.data);
			cudaFree(device_task_cpu.nodes.data);
			cudaFree(device_task_p);

			for(int i=0;i<tmp_netId_to_cellId_s.size();i++){
				cudaFree(tmp_netId_to_cellId_s[i].data);
			}
			cudaFree(device_graph_cpu.netId_to_cellId_s.data);
			for(int i=0;i<G.cellId_to_netId_s.size();i++){
				cudaFree(tmp_cellId_to_netId_s[i].data);
			}
			cudaFree(device_graph_cpu.cellId_to_netId_s.data);
			cudaFree(device_graph_cpu.cellId_to_area.data);
			cudaFree(device_graph_p);
			return move(cost);
		}();
		cost_type best_cost{ G.total_edge_weight,G.total_area,INT_MAX };
		for (int j = 0; j < cost.size(); j++) { //可以用parallel reduction (log(n) complexity)
			if (cost[j] < best_cost) {
				/*cout << "find better move: from cut cost " << best_cost.cut_cost + graph.cut_cost <<
					", balance cost " << best_cost.balance_cost<<" to cut cost " << cost[j].cut_cost+graph.cut_cost
					<< ", balance cost "<< cost[j].balance_cost << endl;*/
				best_cost = cost[j];
			}
		}

		cout << "best cost: cut cost " << best_cost.cut_cost -affected_cut_size+cut_size
			<< ", balance cost " << best_cost.balance_cost
			<< ", idx:" << best_cost.move_idx << endl;
		MOVE(&cellId_to_partId, enumerate_cellId_s, best_cost.move_idx, G);
		auto end = chrono::high_resolution_clock::now();
		auto duration = chrono::duration_cast<chrono::milliseconds>(end - start);
		cout << "round:" << round << ".Time taken: " << duration.count() << " ms" << "\n";
		total_time += duration.count();
	}
	partId_to_partArea =
		cell_id_to_partition_name_generate_partition_name_to_partition_area(cellId_to_partId, G.cellId_to_area);
	cout << "final A area = " << partId_to_partArea[A] << " , B area = " << partId_to_partArea[B] <<
		" ,balance_cost=" << labs(partId_to_partArea[A] - partId_to_partArea[B]) << "\n";
	netId_to_partId_to_cellNum = cell_id_to_partition_name_generate_net_id_to_partition_name_to_cell_num(
		cellId_to_partId, G);
	int cut_size = [&]() {
		int cut_size = 0;
		for (auto& partition_name_to_cell_num_of_net : netId_to_partId_to_cellNum)//O(#nets)
			if (partition_name_to_cell_num_of_net[A] > 0 && partition_name_to_cell_num_of_net[B] > 0)
				cut_size++;
		return cut_size;
		}();
	cout << "final cut size:" << cut_size << "\n";
	cout << "Average time per round: " << total_time / round_num << " ms\n";
	//cost_fout << initial_cut_size << "->" << cut_size << "\n";
	//}
	cost_fout.close();
	cudaStatus = cudaDeviceReset();
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaDeviceReset failed!");
		return 1;
	}
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