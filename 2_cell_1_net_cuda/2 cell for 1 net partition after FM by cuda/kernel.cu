#include "cuda_runtime.h"
#include "device_launch_parameters.h"
using namespace std;
#include <omp.h>
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
        idx = ub = area = lb = 0;
    }
    short idx;
    double  ub , lb;

    unsigned area;
    void report() const {
        printf("part[%d]: area %d, max_area %f, min_area %f\n", idx, area, ub, lb);
    }
};

struct GPart_Edge {
    unsigned idx;
    short weight;
    bool selected = false;
    unsigned node_idx1;
    unsigned node_idx2;
};

struct GPart_Node {
    short idx;
    short area;
    short part_idx;
    short tmp_part_idx = -1;
    std::vector<unsigned> edges;
};

struct GPart_Graph {
    vector<int>cell_program_id_to_input_id;//v[cell_program_id]=cell_input_id
    unordered_map<int, int>cell_input_id_to_program_id;//v[cell_input_id]=cell_program_id

    std::vector<GPart_Partition> parts;
    std::vector<GPart_Node> nodes;
    std::vector<GPart_Edge> edges;
    unsigned total_edge_weight = 0;
    unsigned total_area = 0;
    unsigned cut_cost = 0;

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
    std::vector<unsigned> nodes;
    std::vector<unsigned> edges;
    //std::vector<short> parts;
};



void init_graph(GPart_Graph& graph, const string cell_file, const string net_file, const string partition_file) {

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
        /*for (unsigned i = 0; i < num_nodes; i++) {
            GPart_Node& node = graph.nodes[i];
            node.idx = i;
            node.area = 1 + rand() % max_node_area;
        }*/
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
        //graph.edges.resize(num_edges);    
        /*for (unsigned i = 0; i < num_edges; i++) {
            GPart_Edge& edge = graph.edges[i];
            edge.idx = i;
            edge.weight = 1 + rand() % max_edge_weight;
            edge.node_idx1 = rand() % num_nodes;
            edge.node_idx2 = rand() % num_nodes;
            while (edge.node_idx1 == edge.node_idx2) {
                edge.node_idx2 = rand() % num_nodes;
            }
            graph.nodes[edge.node_idx1].edges.push_back(i);
            graph.nodes[edge.node_idx2].edges.push_back(i);

            graph.total_edge_weight += edge.weight;
        }*/
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
        /*for (unsigned i = 0; i < num_nodes; i++) {
            GPart_Node& node = graph.nodes[i];
            node.part_idx = rand() % num_parts;
            graph.parts[node.part_idx].area += node.area;
        }*/
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
	double lb, ub;
    tie(lb, ub) = []( GPart_Graph& g, const double r, const double epsilon)->pair<int, int> {
        int  max_area = 0;
        g.total_area = 0;
        for (auto& node : g.nodes) {
            g.total_area += node.area;
            if (node.area > max_area) max_area = node.area;
        }
        int  S_opt = ceil(r * g.total_area);//definition of https://chriswalshaw.co.uk/partition/
        double max_S_p = (1 + epsilon) * S_opt;
        return { g.total_area - max_S_p, max_S_p };
        }(graph, 0.5, 0.04);
    for (unsigned i = 0; i < 2; i++) {
        graph.parts[i].idx = i;
        graph.parts[i].lb = lb;
        graph.parts[i].ub = ub;
    }
    for (unsigned i = 0; i < graph.edges.size(); i++) {
        graph.cut_cost += graph.get_cut_cost(i);
    }
}


int random_init_task(GPart_Graph& graph,
    GPart_Task& task, unsigned task_size, ofstream& fout2)
{
    task.nodes.clear();
    task.edges.clear();//report cost
    //task.parts.resize(task_size);
    unsigned num_nodes = graph.nodes.size();

    for (unsigned i = 0; i < task_size; i++) {
        unsigned idx = 0;
        while (1) {
            idx = rand() % num_nodes;
            if (std::find(task.nodes.begin(), task.nodes.end(), idx)
                == task.nodes.end()) {
                task.nodes.push_back(idx);
                break;
            }
        }
    }
	sort(task.nodes.begin(), task.nodes.end());
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
    /*for (int i = 0; i < task_size; i++)
        task.parts[i] = graph.nodes[task.nodes[i]].part_idx;*/
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
/*
double evaluate( GPart_Graph& graph, const GPart_Task& task, unsigned move)
{
    //vector<int> tmp_part_idx(graph.nodes.size(), -1);
    unsigned area[2] = { 0, 0 };
    for (unsigned i = 0; i < task.nodes.size(); i++) {
        GPart_Node& node = graph.nodes[task.nodes[i]];
        //const GPart_Node& node = graph.nodes[task.nodes[i]];
        unsigned part = ((move & (1 << i)) != 0);
        node.tmp_part_idx = part;
        //tmp_part_idx[task.nodes[i]] = part;
        area[part] += node.area;
    }
    double balance_cost = 0;
    for (unsigned i = 0; i < graph.parts.size(); i++) {
        const GPart_Partition& part = graph.parts[i];
        if (part.area + area[i] > part.max_area) {
            //return FLT_MAX;
            //return graph.total_edge_weight;
            return DBL_MAX;//至少有一個移動方法是不動，符合lb,ub之間
        }
        double ratio = double(part.area + area[i]) / double(part.max_area);
        balance_cost += ratio * ratio;//(A area/all area)越遠離0.5,(A area)^2+(all area-A area)^2越大
    }
    unsigned cut_cost = 0;
    for (unsigned edge_idx : task.edges) {
        const GPart_Edge& edge = graph.edges[edge_idx];
        const GPart_Node& n1 = graph.nodes[edge.node_idx1];
        const GPart_Node& n2 = graph.nodes[edge.node_idx2];
        unsigned part1 = (n1.tmp_part_idx != -1) ?n1.tmp_part_idx : n1.part_idx;
        unsigned part2 = (n2.tmp_part_idx != -1) ? n2.tmp_part_idx : n2.part_idx;
        //unsigned part1 = (tmp_part_idx[edge.node_idx1] != -1) ? tmp_part_idx[edge.node_idx1] : n1.part_idx;
        //unsigned part2 = (tmp_part_idx[edge.node_idx2] != -1) ? tmp_part_idx[edge.node_idx2] : n2.part_idx;
        if (part1 != part2) {
            cut_cost += edge.weight;
        }
    }
    return cut_cost + balance_cost / 2;//先比cut size再比平衡,cut size是整數,balance_cost:A area \mapsto ((A area/ub)^2+(B area/ub)^2)/2恆<1
}*/
/*double evaluate(const GPart_Graph& graph, const GPart_Task& task, unsigned move, vector<int>& tmp_part_idx, int task_cut_cost)
{
    unsigned area[2] = { 0, 0 };
    for (unsigned i = 0; i < task.nodes.size(); i++) {
        const GPart_Node& node = graph.nodes[task.nodes[i]];
        unsigned part = ((move & (1 << i)) != 0);
        tmp_part_idx[task.nodes[i]] = part;
        area[part] += node.area;
    }
    double balance_cost = 0;
    for (unsigned i = 0; i < graph.parts.size(); i++) {
        const GPart_Partition& part = graph.parts[i];
        if (part.area + area[i] > part.max_area) {
            //return FLT_MAX;
            return  graph.total_edge_weight;
        }
        double ratio = double(part.area + area[i]) / double(part.max_area);
        balance_cost += ratio * ratio;//(A area/all area)越遠離0.5,(A area)^2+(all area-A area)^2越大
    }
    unsigned cut_cost = 0;
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
    if (cut_cost == task_cut_cost)return  cut_cost + balance_cost / 2;
    else return  cut_cost + balance_cost / 2;

    //return cut_cost + balance_cost/2;//先比cut size再比平衡,cut size是整數,balance_cost:A area \mapsto ((A area/ub)^2+(B area/ub)^2)/2恆<1
}*/
struct device_graph_data{};
struct device_task_data {
    int* edges;
    int nodes_num,edges_num;//可能graph.node.size<max task size
	array<int, MAX_TASK_SIZE> nodes;
	__device__ int find_idx(int node_idx) const{//need sorted array
		int start = 0, end = MAX_TASK_SIZE - 1, mid;
        while (start != end) {
            mid = (start + end) / 2;
            if (nodes[mid] == node_idx)return mid;
            else if (nodes[mid] < node_idx)start = mid + 1;
            else end = mid - 1;
		}
		if (nodes[start] == node_idx)return start;
		return MAX_TASK_SIZE;
};
__device__ int device_evaluate(const device_task_data& task_data, unsigned move) {
    array<int, MAX_TASK_SIZE>& tmp_part_idx;
	tmp_part_idx.fill(-1);
    unsigned area[2] = { 0, 0 };
    for (unsigned i = 0; i < task_data.nodes_num; i++) {
        const GPart_Node& node = graph.nodes[task_data.nodes[i]];
        unsigned part = ((move & (1 << i)) != 0);
        //tmp_part_idx[task.nodes[i]] = part;
		tmp_part_idx[i] = part;
        area[part] += node.area;
    }
    double balance_cost = 0;
    for (unsigned i = 0; i < graph.parts.size(); i++) {
        const GPart_Partition& part = graph.parts[i];
        if (part.area + area[i] > part.max_area) {
            //return FLT_MAX;
            return  graph.total_edge_weight;
        }
        double ratio = double(part.area + area[i]) / double(part.max_area);
        balance_cost += ratio * ratio;//(A area/all area)越遠離0.5,(A area)^2+(all area-A area)^2越大
    }
    unsigned cut_cost = 0;
    for (unsigned edge_idx : task_data.edges) {
        const GPart_Edge& edge = graph.edges[edge_idx];
        const GPart_Node& n1 = graph.nodes[edge.node_idx1];
        const GPart_Node& n2 = graph.nodes[edge.node_idx2];
		int idx;
        //unsigned part1 = (tmp_part_idx[edge.node_idx1] != -1) ? tmp_part_idx[edge.node_idx1] : n1.part_idx;
        unsigned part1 = ((idx= task_data.find_idx(edge.node_idx1))<MAX_TASK_SIZE) ? tmp_part_idx[idx] : n1.part_idx;
        //unsigned part2 = (tmp_part_idx[edge.node_idx2] != -1) ? tmp_part_idx[edge.node_idx2] : n2.part_idx;
        unsigned part2 = ((idx = task_data.find_idx(edge.node_idx2)) < MAX_TASK_SIZE) ? tmp_part_idx[idx] : n2.part_idx;
        if (part1 != part2) {
            cut_cost += edge.weight;
        }
    }
    return  cut_cost + balance_cost / 2;

}
__global__  void evaluate_kernel( /* device graph data */, /* device task data */,
    int task_cut_cost, int task_size, double* cost, int num_explores ) {
    int j = blockIdx.x * blockDim.x + threadIdx.x;

    // 這裡需將 evaluate 的邏輯搬到 device 上
    // 你需要將 graph, task, tmp_part_idx 等資料結構轉成 device-friendly 格式
    // 這裡僅示意
	unsigned stride = blockDim.x * gridDim.x;
    for (j = 0; j < num_explores; j += stride) {
        cost[j] = device_evaluate(/* device graph, device task, j, ... */);
    }
}
void commit(GPart_Graph& graph, GPart_Task& task, unsigned best_move)
{
    for (unsigned i = 0; i < task.nodes.size(); i++) {
        GPart_Node& node = graph.nodes[task.nodes[i]];
        node.part_idx = ((best_move & (1 << i)) != 0);
        node.tmp_part_idx = -1;
    }

    for (unsigned idx : task.edges) {
        graph.edges[idx].selected = false;
    }
}

int main(int argc, char* argv[]) {

    /*
    if(argc <= 2 || !strcmp(argv[1], "-help")) {
      help_command();
      return 0;
    }
    */
    GPart_Graph graph;
    string file_name_base =/*"../../simple_testcase/simple_graph_testcase";//*/"../../testcase/delaunay_n10";
    init_graph(graph, file_name_base + ".cells", file_name_base + ".nets", file_name_base + ".partitions");
    ofstream fout2("log.txt");

    printf("Initi partition:\n");
    profile_graph(graph);

    unsigned num_trails = 5;
    int task_size = min(MAX_TASK_SIZE, (int)graph.nodes.size()); //smaller than 32
    const int round_num = 10;double total_time = 0;
    GPart_Task task;
    for (unsigned round = 1; round <= num_trails; round++) {
        auto wcts = std::chrono::system_clock::now();
        int task_cut_cost = random_init_task(graph, task, task_size, fout2);
        update_graph(graph, task, false);

        unsigned num_explores = 1 << task_size;
        //printf("num explores %d\n", num_explores);
        vector<int>tmp_part_idx(graph.nodes.size(), -1);
/*#pragma omp parallel num_threads(4) firstprivate(local_min,local_min_idx,tmp_part_idx,local_out_of_lb_ub_num,local_cost_not_changed_num)
        {
#pragma omp for schedule(dynamic, 10000)
            for (int j = 0; j < num_explores; j++) {
                cost[j] = evaluate(graph, task, j, tmp_part_idx, task_cut_cost);
                for (unsigned i = 0; i < task.nodes.size(); i++) {
                    tmp_part_idx[task.nodes[i]] = -1;
                }

                if (cost < local_min) {
                    local_min = cost;
                    local_min_idx = j;
                }
            }
#pragma omp critical
            {
                if (local_min < best_cost) {
                    best_cost = local_min;
                    best_move = local_min_idx;
                }
            }
        }*/

        int threadsPerBlock = 512;
        int blocksPerGrid = 16;
        double* d_cost;
        cudaMalloc(&d_cost, num_explores * sizeof(double));

        // 準備 device graph, device task, 以及其他必要資料
        // 你需要將 graph, task 等資料複製到 device（建議用 struct of arrays 或 flat array）

        evaluate_kernel << <blocksPerGrid, threadsPerBlock >> > (
            /* device graph data */,
            /* device task data */,
            task_cut_cost,
            task_size,
            d_cost,
            num_explores
            );
        cudaDeviceSynchronize();

        // 將結果複製回 host
        double* cost = new double[num_explores];
        cudaMemcpy(cost, d_cost, num_explores * sizeof(double), cudaMemcpyDeviceToHost);

        // 找出 best_move
        int best_move = -1;
        double best_cost = DBL_MAX;
        for (int j = 0; j < num_explores; j++) {
            //parallel reduction (if thread infinity it has log(n) complexity)
            if (cost[j] < best_cost) {
                best_cost = cost[j];
                best_move = j;
            }
        }

        cudaFree(d_cost);
		delete[] cost;
        fout2 << "the rate of move not change cost: " << gobal_cost_not_changed_num * 1.0 / num_explores << endl;
        fout2 << "the rate of move out of lb/ub: " << gobal_out_of_lb_ub_num * 1.0 / num_explores << endl;
        commit(graph, task, best_move);

        update_graph(graph, task, true);

        std::chrono::duration<double> wctduration =
            (std::chrono::system_clock::now() - wcts);
        printf("\nRound %d (%.2lf seconds)\n", round, wctduration.count());
        total_time += wctduration.count();
        profile_graph(graph);
        fflush(stdout);
    }
    fout2.close();
    printf("average time: %.2lf seconds\n", total_time / round_num);
    return 1;
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