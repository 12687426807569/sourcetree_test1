#include "cuda_runtime.h"
#include "device_launch_parameters.h"
using namespace std;
//#include <omp.h>
#include <vector>
#include<tuple>
#include <algorithm>
#include <iostream>
#include <cfloat>
#include <chrono>
#include<fstream>
#include<array>
#include<cmath>
#include<unordered_map>
#include<string>
constexpr int A = 0, B = 1;
const int MAX_TASK_SIZE = 20;
struct GPart_Partition {
    GPart_Partition() {
        idx = ub = area =  0;
    }
    short idx;
    unsigned  ub ;//double不准

    unsigned area;
    void report() const {
        printf("part[%d]: area %d, max_area %d\n", idx, area, ub);
    }
};

struct GPart_Edge {
    unsigned idx;
    unsigned weight;
    bool selected = false;
    unsigned node_idx1;
    unsigned node_idx2;
};

struct GPart_Node {
    short idx;
    short area;
    short part_idx;
    //short tmp_part_idx = -1;
    std::vector<unsigned> edges;
};

struct GPart_Graph {
    vector<int>cell_program_id_to_input_id;//a[cell_program_id]=cell_input_id
    unordered_map<int, int>cell_input_id_to_program_id;//a[cell_input_id]=cell_program_id

    std::vector<GPart_Partition> parts;
    std::vector<GPart_Node> nodes;
    std::vector<GPart_Edge> edges;
    int total_edge_weight = 0;
    int total_area = 0;
    int cut_cost = 0;

    unsigned get_cut_cost(unsigned edge_idx) const {
        const GPart_Edge& edge = edges[edge_idx];
        if (nodes[edge.node_idx1].part_idx
            != nodes[edge.node_idx2].part_idx) {
            return edge.weight;
        }
        return 0;
    }
};

struct GPart_Task {
    std::vector<int> nodes;
    std::vector<unsigned> edges;
};

void init_graph(GPart_Graph& graph, const string&cell_file, const string&net_file, const string&partition_file) {
    const unsigned num_parts = 2;
    graph.parts.resize(num_parts);
    GPart_Graph& G = graph;
    {
        ifstream fin;
        fin.open(cell_file);
        string cell_id, area;
        unordered_map<int, int>cell_input_id_to_area;//v[cell_input_id]=area
        while (fin >> cell_id)//O(#cells)
        {
            fin >> area;
            cell_input_id_to_area[stoi(cell_id.substr(1))] = stoi(area);//O(1))
        }
        graph.nodes.resize(cell_input_id_to_area.size());
        int i = 0;
        for (auto pair_element : cell_input_id_to_area) {//O(#cells)
            int cell_input_id= pair_element.first, area= pair_element.second;
            GPart_Node& node = graph.nodes[i];
            node.idx = i;
            node.area = area;
            G.cell_program_id_to_input_id.push_back(cell_input_id);
            G.cell_input_id_to_program_id[cell_input_id] = i;//O(1)
            i++;
        }
        fin.close();

        fin.open(net_file);
        string t;
        i = 0;
        while (fin >> t)//O(maxDegree * #nets)
        {
            graph.edges.push_back(GPart_Edge());
            GPart_Edge& edge = graph.edges[i];
            edge.idx = i;
            edge.weight = 1;
            fin >> t >> t;
            vector<int> cells;
            while (fin >> t)//O(maxDegree)
            {
                if (t == "}")break;
                int cell_id = G.cell_input_id_to_program_id[stoi(t.substr(1))];
                cells.push_back(cell_id); //O(1)
            }
            edge.node_idx1 = cells[0];
            edge.node_idx2 = cells[1];
            graph.nodes[edge.node_idx1].edges.push_back(i);
            graph.nodes[edge.node_idx2].edges.push_back(i);
            graph.total_edge_weight += edge.weight;

            i++;
        }
        fin.close();
    }
    {
        string t;
        ifstream fin;
        fin.open(partition_file);
        fin >> t >> t;
        int cut_size = stoi(t);
        fin >> t >> t;
        for (int num = stoi(t); num > 0; num--) {//O(#cells)
            fin >> t;
            int cell_id = G.cell_input_id_to_program_id.at(stoi(t.substr(1)));//O(1)
            GPart_Node& node = graph.nodes[cell_id];
            node.part_idx = A;
            graph.parts[node.part_idx].area += node.area;
        }
        fin >> t >> t;
        for (int num = stoi(t); num > 0; num--) {//O(#cells)
            fin >> t;
            int cell_id = G.cell_input_id_to_program_id.at(stoi(t.substr(1)));//O(1)
            GPart_Node& node = graph.nodes[cell_id];
            node.part_idx = B;
            graph.parts[node.part_idx].area += node.area;
        }
        fin.close();
    }
	int lb, ub;
    tie(lb, ub) = [&graph]( const /*double r*/unsigned P, const /*double epsilon*/unsigned x)->pair<int, int> {
        int  max_area = 0;
        graph.total_area = 0;
        for (auto& node : graph.nodes) {
            graph.total_area += node.area;
            if (node.area > max_area) max_area = node.area;
        }
        //definition of https://chriswalshaw.co.uk/partition/
        //int  S_opt = ceil(r * graph.total_area);//double不准
        int  S_opt = (graph.total_area + P - 1) / P;
        //int max_S_p =floor( (1 + epsilon) * S_opt);
        int max_S_p = S_opt * (100 + x) / 100;
        return { graph.total_area - max_S_p, max_S_p };
        }(/* 0.5*/2, /*0.1*/10);
    for (unsigned i = 0; i < 2; i++) {
        graph.parts[i].idx = i;
        graph.parts[i].ub = ub;
    }
    for (unsigned i = 0; i < graph.edges.size(); i++) {
        graph.cut_cost += graph.get_cut_cost(i);
    }
}

int random_init_task(GPart_Graph& graph,
    GPart_Task& task, unsigned task_size, ofstream& fout2){
    task.nodes.clear();
    task.edges.clear();//report cost
    //task.parts.resize(task_size);
    unsigned num_nodes = graph.nodes.size();

    for (unsigned i = 0; i < task_size; i++) {
        unsigned idx = 0;
        while (1) {
            idx = rand() % num_nodes;
			if (std::find(task.nodes.begin(), task.nodes.end(), idx)//可以用unordered_set,find in O(1)，
                == task.nodes.end()) {
                task.nodes.push_back(idx);
                break;
            }
        }
    }
	//sort(task.nodes.begin(), task.nodes.end());//find idx in O(logn) by binary search
    {
        struct cout_part_format {
            int cell_num = 0;
            int area = 0;
        }a[2];
        for (auto cell : task.nodes) {
            a[graph.nodes[cell].part_idx].cell_num++;
            a[graph.nodes[cell].part_idx].area += graph.nodes[cell].area;
        }
        fout2 << "this task part A: cell num " << a[0].cell_num << ", area " << a[0].area
            << "; part B: cell num " << a[1].cell_num << ", area " << a[1].area << endl;
    }
    for (unsigned i = 0; i < task_size; i++) {
        for (unsigned idx : graph.nodes[task.nodes[i]].edges) {
            if (!graph.edges[idx].selected) {
                task.edges.push_back(idx);
                graph.edges[idx].selected = true;
            }
        }
    }
    int task_cut_cost = 0;
    for (auto edge : task.edges) {
        task_cut_cost += graph.get_cut_cost(edge);
    }
    fout2 << "this task cut cost " << task_cut_cost << endl;
    return task_cut_cost;
}

void update_graph(GPart_Graph& graph, const GPart_Task& task, bool add)
{
    for (unsigned idx : task.nodes) {
        const GPart_Node& node = graph.nodes[idx];
        GPart_Partition& part = graph.parts[node.part_idx];
        if (add) {
            part.area += node.area;
        }
        else {
            part.area -= node.area;
        }
    }

    //printf("edge size = %ld\n", task.edges.size());
    for (unsigned idx : task.edges) {
        //if(graph.get_cut_cost(idx)>0) cout << "edge idx " << idx << " cut cost " << graph.get_cut_cost(idx) << endl;
        if (add) {
            graph.cut_cost += graph.get_cut_cost(idx);
        }
        else {
            graph.cut_cost -= graph.get_cut_cost(idx);
        }
    }
}

void profile_graph(const GPart_Graph& graph)
{
    for (unsigned i = 0; i < graph.parts.size(); i++) {
        graph.parts[i].report();
    }

    printf("Total Edge Weight %d\n", graph.total_edge_weight);
    printf("Total Cut Cost %d\n", graph.cut_cost);
}

struct cost_type {
    int cut_cost;
    int balance_cost;
    int move_idx;
    __host__ __device__ bool operator<(const cost_type&right) {
        if (cut_cost < right.cut_cost)return true;
        if (cut_cost == right.cut_cost && balance_cost < right.balance_cost)return true;
		return false;
    }
};
cost_type evaluate(const GPart_Graph& graph, const GPart_Task& task, int move_idx, vector<int>& tmp_part_idx, 
    int task_cut_cost){
    unsigned area[2] = { 0, 0 };
    for (unsigned i = 0; i < task.nodes.size(); i++) {
        const GPart_Node& node = graph.nodes[task.nodes[i]];
        unsigned part = ((move_idx & (1 << i)) != 0);
        tmp_part_idx[task.nodes[i]] = part;
        area[part] += node.area;
    }
    for (unsigned i = 0; i < graph.parts.size(); i++) {
        const GPart_Partition& part = graph.parts[i];
        if (part.area + area[i] > part.ub) {
            return  { graph.total_edge_weight ,graph.total_area,move_idx };
        }
    }
    int balance_cost = labs(graph.parts[A].area + area[A] - graph.parts[B].area - area[B]);
    int cut_cost = 0;
    for (unsigned edge_idx : task.edges) {
        const GPart_Edge& edge = graph.edges[edge_idx];
        const GPart_Node& n1 = graph.nodes[edge.node_idx1];
        const GPart_Node& n2 = graph.nodes[edge.node_idx2];
        //unsigned part1 = (n1.tmp_part_idx != -1) ?n1.tmp_part_idx : n1.part_idx;
        //unsigned part2 = (n2.tmp_part_idx != -1) ? n2.tmp_part_idx : n2.part_idx;
        unsigned part1 = (tmp_part_idx[edge.node_idx1] != -1) ? tmp_part_idx[edge.node_idx1] : n1.part_idx;
        unsigned part2 = (tmp_part_idx[edge.node_idx2] != -1) ? tmp_part_idx[edge.node_idx2] : n2.part_idx;
        if (part1 != part2) {
            cut_cost += edge.weight;
        }
    }
    return  {cut_cost, balance_cost,move_idx };
}
void evaluate_array(vector<cost_type>& cost, const GPart_Graph& graph, const GPart_Task& task, 
	int task_cut_cost, unsigned num_explores) {//openmp parallel for
    vector<int>tmp_part_idx(graph.nodes.size(), -1);
    int cpu_thread_num = 1, thread_idx = 0, move_idx = thread_idx;
    cost[thread_idx] = { INT_MAX,INT_MAX,0 };
    for (; move_idx < num_explores; move_idx += cpu_thread_num) {
        auto now_cost =  evaluate(graph, task, move_idx, tmp_part_idx, task_cut_cost);
        for (unsigned i = 0; i < task.nodes.size(); i++) {
            tmp_part_idx[task.nodes[i]] = -1;
        }
        if (now_cost<cost[thread_idx]  )cost[thread_idx] = now_cost;
    }
}
struct device_graph{
    struct device_graph_node {
        int area, part_idx;
    }*nodes;
    struct device_graph_edge {
        unsigned node_idx1, node_idx2,weight;
    }*edges;
    struct device_graph_part {
        unsigned area;
        unsigned ub;
    } parts[2];
    int total_edge_weight, total_area;
};
struct device_task { 
    int*nodes,* edges;
    int nodes_num, edges_num;
    struct hash_node {
        int node_idx, task_idx, next_idx;
    }*node_idx_to_task_idx_data;//inv_nodes(nodes[i])=i,inv_nodes(not in nodes)=nodes_num
    /*__device__ int find_idx(int node_idx) const {//need sorted array
        int start = 0, end = nodes_num - 1, mid;
        while (start != end) {
            mid = (start + end) / 2;
            if (nodes[mid] == node_idx)return mid;
            else if (nodes[mid] < node_idx)start = mid + 1;
            else end = mid - 1;
        }
        if (nodes[start] == node_idx)return start;
        return nodes_num;
    }*/
    __device__ int node_idx_to_task_idx(int node_idx) const {
		int hash_idx = node_idx % (2*nodes_num);
		//linked list,data[0~2*nodes_num-1} are dummy heads,data[2*nodes_num~] are list nodes,load factor=1/2
		hash_idx = node_idx_to_task_idx_data[hash_idx].next_idx;
        while (hash_idx !=-1) {
			if (node_idx_to_task_idx_data[hash_idx].node_idx == node_idx)return node_idx_to_task_idx_data[hash_idx].task_idx;
            hash_idx = node_idx_to_task_idx_data[hash_idx].next_idx;
        }
		return nodes_num;
	}
};
__device__ cost_type device_evaluate(const device_graph& graph, const device_task& task, int move_idx,
    int tmp_part_idx[], int task_cut_cost) {
    unsigned area[2] = { 0, 0 };
    for (unsigned i = 0; i < task.nodes_num; i++) {
        const auto& node = graph.nodes[task.nodes[i]];
        unsigned part = ((move_idx & (1 << i)) != 0);
        //tmp_part_idx[task.nodes[i]] = part;
        tmp_part_idx[i] = part;
        area[part] += node.area;
    }
    for (unsigned i = 0; i < 2; i++) {
        const auto& part = graph.parts[i];
        if (part.area + area[i] > part.ub) {
            return  { graph.total_edge_weight ,graph.total_area ,move_idx };
        }
    }
    int balance_cost = labs(graph.parts[A].area + area[A] - graph.parts[B].area - area[B]);
    int cut_cost = 0;
    for (int i = 0; i < task.edges_num; i++) {
        unsigned edge_idx = task.edges[i];
        const auto& edge = graph.edges[edge_idx];
        const auto& n1 = graph.nodes[edge.node_idx1];
        const auto& n2 = graph.nodes[edge.node_idx2];
        int idx;
        //unsigned part1 = (tmp_part_idx[edge.node_idx1] != -1) ? tmp_part_idx[edge.node_idx1] : n1.part_idx;
        unsigned part1 = ((idx = task.node_idx_to_task_idx(edge.node_idx1)) < task.nodes_num) ? tmp_part_idx[idx] : n1.part_idx;
        //unsigned part2 = (tmp_part_idx[edge.node_idx2] != -1) ? tmp_part_idx[edge.node_idx2] : n2.part_idx;
        unsigned part2 = ((idx = task.node_idx_to_task_idx(edge.node_idx2)) < task.nodes_num) ? tmp_part_idx[idx] : n2.part_idx;
        if (part1 != part2) 
            cut_cost += edge.weight;
    }
    return  { cut_cost, balance_cost,move_idx };
}
__global__  void evaluate_array_kernel(cost_type* cost, device_graph* device_graph1_p, device_task* device_task1_p
    , int task_cut_cost, int num_explores) {
    unsigned cost_idx = blockIdx.x * blockDim.x + threadIdx.x,  move_idx = cost_idx;

    unsigned stride = blockDim.x * gridDim.x;
	int tmp_part_idx[MAX_TASK_SIZE];
    cost[cost_idx] = { INT_MAX,INT_MAX,0 };
    for (; move_idx < num_explores; move_idx += stride) {
        for (int i = 0; i < device_task1_p->nodes_num; i++)tmp_part_idx[i] = -1;
        cost_type now_cost = device_evaluate(*device_graph1_p, *device_task1_p, move_idx, tmp_part_idx, task_cut_cost);
        if ( now_cost<cost[cost_idx] )cost[cost_idx] = now_cost;
    }
}
void commit(GPart_Graph& graph, GPart_Task& task, unsigned best_move)
{
    for (unsigned i = 0; i < task.nodes.size(); i++) {
        GPart_Node& node = graph.nodes[task.nodes[i]];
        node.part_idx = ((best_move & (1 << i)) != 0);
        //node.tmp_part_idx = -1;
    }

    for (unsigned idx : task.edges) {
        graph.edges[idx].selected = false;
    }
}

int main(int argc, char* argv[]) {
    bool use_cuda =1;
    cudaError_t cudaStatus = cudaSetDevice(0);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaSetDevice failed!  Do you have a CUDA-capable GPU installed?");
		exit(1);
    }

    GPart_Graph graph;
    string file_name_base =/*"../../simple_testcase/simple_graph_testcase";//*/"../../testcase/delaunay_n10";
    init_graph(graph, file_name_base + ".cells", file_name_base + ".nets", file_name_base + ".partitions");
    ofstream fout2("log.txt");

    printf("Initi partition:\n");
    profile_graph(graph);

    unsigned num_trails =10;
    int task_size = min(MAX_TASK_SIZE, (int)graph.nodes.size());
    double total_time = 0;
    GPart_Task task;
    for (unsigned round = 1; round <= num_trails; round++) {
        auto wcts = std::chrono::system_clock::now();
        int task_cut_cost = random_init_task(graph, task, task_size, fout2);
        update_graph(graph, task, false);

        unsigned num_explores = 1 << task_size;
        //printf("num explores %d\n", num_explores);
        vector<cost_type> cost;
        //evaluate_cost_array
        //[&cost, &graph, &task, &tmp_part_idx, task_cut_cost, num_explores, use_cuda]()
        {
            if (use_cuda) {
                int threadsPerBlock = 512;
                int blocksPerGrid = 16;
				cost.resize(blocksPerGrid * threadsPerBlock);
                cost_type* dev_cost;
                cudaMalloc(&dev_cost, blocksPerGrid * threadsPerBlock * sizeof(cost_type));

                device_task*device_task1_p;
				cudaMalloc(&device_task1_p, sizeof(device_task));
				device_task device_task_cpu;

                cudaMalloc(&device_task_cpu.nodes, task.nodes.size() * sizeof(int));
                cudaMemcpy(device_task_cpu.nodes, task.nodes.data(), task.nodes.size() * sizeof(int), cudaMemcpyHostToDevice);
                device_task_cpu.nodes_num = task.nodes.size();
				
                device_task::hash_node* hash_data = new device_task::hash_node[3 * task.nodes.size()];
				for (int i = 0; i <  task.nodes.size(); i++)hash_data[i] = device_task::hash_node{ -1,-1,-1 };
                for (int i = 0, free_idx = 2 * task.nodes.size();i < task.nodes.size(); i++) {
                    int hash_idx = task.nodes[i] % (2 * task.nodes.size());
					//list head is hash_data[0~2*task.nodes.size()-1], list node is hash_data[2*task.nodes.size()~],load factor=1/2
                    hash_data[free_idx] = device_task::hash_node{ task.nodes[i], i,hash_data[hash_idx].next_idx };
                    hash_data[hash_idx].next_idx = free_idx;
                    free_idx++;
                }
				cudaMalloc(&device_task_cpu.node_idx_to_task_idx_data, 3* task.nodes.size() * sizeof(device_task::hash_node));
				cudaMemcpy(device_task_cpu.node_idx_to_task_idx_data, hash_data, 3 * task.nodes.size() * sizeof(device_task::hash_node), cudaMemcpyHostToDevice);
                delete[] hash_data;

                cudaMalloc(&device_task_cpu.edges, task.edges.size() * sizeof(int));
                cudaMemcpy(device_task_cpu.edges, task.edges.data(), task.edges.size() * sizeof(int), cudaMemcpyHostToDevice);
                device_task_cpu.edges_num = task.edges.size();

                cudaMemcpy(device_task1_p, &device_task_cpu, sizeof(device_task), cudaMemcpyHostToDevice);

                device_graph* device_graph1_p;
				cudaMalloc(&device_graph1_p, sizeof(device_graph));
				device_graph device_graph_cpu;

				device_graph::device_graph_node* device_graph_cpu_nodes = new device_graph::device_graph_node[graph.nodes.size()];
                for (int i = 0; i < graph.nodes.size(); i++) 
                    device_graph_cpu_nodes[i] = device_graph::device_graph_node{ graph.nodes[i].area,graph.nodes[i].part_idx };				
				cudaMalloc(&device_graph_cpu.nodes, graph.nodes.size() * sizeof(device_graph::device_graph_node));
                cudaMemcpy(device_graph_cpu.nodes, device_graph_cpu_nodes, graph.nodes.size() * sizeof(device_graph::device_graph_node), cudaMemcpyHostToDevice);
                delete[] device_graph_cpu_nodes;
                device_graph_cpu.total_area = graph.total_area;
                device_graph_cpu.total_edge_weight = graph.total_edge_weight;
                for (int i = 0; i < 2; i++) {
                    device_graph_cpu.parts[i] = device_graph::device_graph_part{ graph.parts[i].area,graph.parts[i].ub };
				}

				device_graph::device_graph_edge* device_graph_cpu_edges = new device_graph::device_graph_edge[graph.edges.size()];
                for (int i = 0; i < graph.edges.size(); i++)
					device_graph_cpu_edges[i] = device_graph::device_graph_edge{ graph.edges[i].node_idx1,graph.edges[i].node_idx2,graph.edges[i].weight };
                cudaMalloc(&device_graph_cpu.edges, graph.edges.size() * sizeof(device_graph::device_graph_edge));
                cudaMemcpy(device_graph_cpu.edges, device_graph_cpu_edges, graph.edges.size() * sizeof(device_graph::device_graph_edge), cudaMemcpyHostToDevice);
                delete[] device_graph_cpu_edges;
                cudaMemcpy(device_graph1_p, &device_graph_cpu, sizeof(device_graph), cudaMemcpyHostToDevice);
				
                evaluate_array_kernel << <blocksPerGrid, threadsPerBlock >> > (
                    dev_cost,
                    device_graph1_p,
                    device_task1_p,
                    task_cut_cost,
					num_explores
                    );
                cudaDeviceSynchronize();
				cudaMemcpy(cost.data(), dev_cost, blocksPerGrid * threadsPerBlock * sizeof(cost_type), cudaMemcpyDeviceToHost);
                cudaFree(dev_cost);

				cudaFree(device_task_cpu.node_idx_to_task_idx_data);
				cudaFree(device_task_cpu.edges);
				cudaFree(device_task_cpu.nodes);
				cudaFree(device_task1_p);

				cudaFree(device_graph_cpu.edges);
				cudaFree(device_graph_cpu.nodes);
				cudaFree(device_graph1_p);

            }
            else {
                int cpu_thread_num = 1;
				cost.resize(cpu_thread_num);
                evaluate_array(cost, graph, task, task_cut_cost, num_explores);
            }
        }

        // 找出 best_move
        cost_type best_cost{ INT_MAX,INT_MAX ,0};
        for (int j = 0; j < cost.size(); j++) { //可以用parallel reduction (log(n) complexity)
			if (cost[j] < best_cost) {
				cout << "find better move: cut cost " << cost[j].cut_cost << ", balance cost " << cost[j].balance_cost << endl;
                best_cost = cost[j];
            }
        }
        cout<< "best cost: cut cost " << best_cost.cut_cost
			<< ", balance cost " << best_cost.balance_cost << endl;

        commit(graph, task, best_cost.move_idx);

        update_graph(graph, task, true);

        std::chrono::duration<double> wctduration =
            (std::chrono::system_clock::now() - wcts);
        printf("\nRound %d (%.2lf seconds)\n", round, wctduration.count());
        total_time += wctduration.count();
        profile_graph(graph);
        fflush(stdout);
    }
    fout2.close();
    printf("average time: %.2lf seconds\n", total_time / num_trails);
    cudaStatus = cudaDeviceReset();
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaDeviceReset failed!");
        return 1;
    }
    return 0;
}/*
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
}*/