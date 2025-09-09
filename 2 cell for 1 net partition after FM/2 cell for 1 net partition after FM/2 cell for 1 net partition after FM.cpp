using namespace std;
#include <omp.h>
#include <vector>
#include <algorithm>
#include <iostream>
#include <cfloat>
#include <chrono>
#include<fstream>
#include<cmath>
#include<tuple>
#include<unordered_map>
#include<string>
#include<numeric>
#include<random>
constexpr int A = 0, B = 1;
struct GPart_Partition {
    GPart_Partition() {
        idx  = area = lb = 0;
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
    vector<int>cell_program_id_to_input_id;//a[cell_program_id]=cell_input_id
    unordered_map<int, int>cell_input_id_to_program_id;//a[cell_input_id]=cell_program_id

    std::vector<GPart_Partition> parts;
    std::vector<GPart_Node> nodes;
    std::vector<GPart_Edge> edges;
    unsigned total_edge_weight = 0;
    unsigned total_area = 0;
    unsigned cut_cost = 0;

    unsigned get_cut_cost(unsigned edge_idx) const {
        const GPart_Edge& edge = edges[edge_idx];
        if (nodes[edge.node_idx1].part_idx != nodes[edge.node_idx2].part_idx) {
            return edge.weight;
        }
        return 0;
    }
};
struct GPart_Task {
    std::vector<unsigned> nodes;
    std::vector<unsigned> edges;
    std::vector<short> parts;
};
void init_graph(GPart_Graph& graph, const string cell_file, const string net_file, const string partition_file) {

    const unsigned num_parts = 2;

    graph.parts.resize(num_parts);
    GPart_Graph& G = graph;
    {
        ifstream fin;
        fin.open(cell_file);
        string cell_id, area;
        unordered_map<int, int>cell_input_id_to_area;//a[cell_input_id]=area
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
        for (auto [cell_input_id, area] : cell_input_id_to_area) {//O(#cells)
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
    const auto [lb, ub] = [](GPart_Graph&g, const double r, const double epsilon)->pair<int, int> {
        int   max_area = 0;
        g.total_area = 0;
        for (auto& node : g.nodes) {
            g.total_area += node.area;
            if (node.area > max_area) max_area = node.area;
        }
        int  S_opt = ceil(r * g.total_area);//definition of https://chriswalshaw.co.uk/partition/
        int max_S_p =floor( (1 + epsilon) * S_opt);
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
int random_init_task(GPart_Graph& graph, GPart_Task& task, int&task_size, ofstream& fout2/*, vector<int>& range*/)
{
    task.nodes.clear();
    task.edges.clear();

	//task_size = min(task_size, (int)range.size());
	//task.nodes.insert(task.nodes.end(),make_move_iterator( range.end()- task_size), make_move_iterator(range.end()) );
	//range.erase(range.end() - task_size, range.end());

    task.parts.resize(task_size);
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
    for(int i=0;i<task_size;i++)
		task.parts[i] = graph.nodes[task.nodes[i]].part_idx;
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
tuple<int,int> evaluate(const GPart_Graph& graph, const GPart_Task& task, unsigned move,
    vector<int>& tmp_part_idx){
	//use int cost instead of double cost to avoid precision problem
    unsigned area[2] = { 0, 0 };
    for (unsigned i = 0; i < task.nodes.size(); i++) {
        //GPart_Node& node = graph.nodes[task.nodes[i]];
        const GPart_Node& node = graph.nodes[task.nodes[i]];
        unsigned part = ((move & (1 << i)) != 0);
        //node.tmp_part_idx = part;
        tmp_part_idx[task.nodes[i]] = part;
        area[part] += node.area;
    }
    for (unsigned i = 0; i < graph.parts.size(); i++) {
        const GPart_Partition& part = graph.parts[i];
        if (part.area + area[i] > part.ub) {
            //return FLT_MAX;
            return { graph.total_edge_weight,graph.total_area};
        }
    }
	int imbalance_cost = labs(graph.parts[A].area + area[A] - graph.parts[B].area - area[B]);
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
    return { cut_cost ,imbalance_cost  };

    //return cut_cost + balance_cost/2;//先比cut size再比平衡,cut size是整數,balance_cost:A area \mapsto ((A area/ub)^2+(B area/ub)^2)/2恆<1
}
void commit(GPart_Graph& graph, GPart_Task& task, unsigned best_move)
{
    for (unsigned i = 0; i < task.nodes.size(); i++) {
        GPart_Node& node = graph.nodes[task.nodes[i]];
        if (node.part_idx != ((best_move & (1 << i)) != 0)) cout << "node idx " << node.idx << " from part " << node.part_idx << " to part " << ((best_move & (1 << i)) != 0) << endl;
        node.part_idx = ((best_move & (1 << i)) != 0);
        node.tmp_part_idx = -1;
    }

    for (unsigned idx : task.edges) {
        graph.edges[idx].selected = false;
    }
}
int main(int argc, char* argv[]) {
    GPart_Graph graph;
    string file_name_base = /*"../../simple_testcase/simple_graph_testcase";//*/"../../testcase/delaunay_n10";
    init_graph(graph, file_name_base + ".cells", file_name_base + ".nets", file_name_base + ".partitions");
    ofstream fout2("log.txt");

    printf("Initi partition:\n");
    profile_graph(graph);

    constexpr unsigned num_trails = 1;
    int task_size = min(20, (int)graph.nodes.size()); 
    double total_time = 0;
    GPart_Task task;
    //vector<int> range(graph.nodes.size());
    //iota(range.begin(), range.end(), 0);
    //shuffle(range.begin(), range.end(), mt19937(12345));
    for (unsigned round = 1; round <= num_trails; round++) {
        auto wcts = std::chrono::system_clock::now();
        int task_cut_cost = random_init_task(graph, task, task_size, fout2/*, range*/);
        /*int not_changed_idx = 0;
		for (unsigned i = 0; i < task_size; i++)
			not_changed_idx += task.parts[i] << (i);*/
        update_graph(graph, task, false);

        unsigned num_explores = 1 << task_size;
		ofstream cut_distribution("cut_distribution.txt"),imbalance_distribution("imbalance_distribution.txt");
        //printf("num explores %d\n", num_explores);
        int best_move = -1, local_min_idx = -1, local_out_of_lb_ub_num = 0, gobal_out_of_lb_ub_num = 0,
            local_cost_not_changed_num = 0, gobal_cost_not_changed_num = 0;
        int best_cut_cost = INT_MAX, best_imbalance_cost = INT_MAX, local_min_cut_cost = INT_MAX,local_min_imbalance_cost=INT_MAX;
        vector<int>tmp_part_idx(graph.nodes.size(), -1);
//#pragma omp parallel num_threads(4) firstprivate(local_min,local_min_idx,tmp_part_idx,local_out_of_lb_ub_num,local_cost_not_changed_num)
        {
//#pragma omp for schedule(dynamic, 10000)
            for (int j = 0; j < num_explores; j++) {
                /*if (j == not_changed_idx)
                    cout << "test_not_change\n";*/
                auto [cut_cost, imbalance_cost] = evaluate(graph, task, j, tmp_part_idx);
				/*cut_distribution << cut_cost << endl;
				imbalance_distribution << imbalance_cost << endl;*/
                for (unsigned i = 0; i < task.nodes.size(); i++) {//cuda版本tmp_part_idx大小可設task.nodes.size()來節省gputhread記憶體
                    tmp_part_idx[task.nodes[i]] = -1;
                }
                /*if (imbalance_cost == graph.total_area) {
                    local_out_of_lb_ub_num++;
                    continue;
                }*/

                if (cut_cost < local_min_cut_cost) {
                    local_min_cut_cost = cut_cost;
                    local_min_imbalance_cost = imbalance_cost;
                    local_min_idx = j;
                }
                else if (cut_cost == local_min_cut_cost && imbalance_cost < local_min_imbalance_cost) {
                    local_min_imbalance_cost = imbalance_cost;
                    local_min_idx = j;
                }
                /*if (task_cut_cost == cut_cost)
                    local_cost_not_changed_num++;*/
            }
//#pragma omp critical
            {
                if (local_min_cut_cost < best_cut_cost) {
                    best_cut_cost = local_min_cut_cost;
					best_imbalance_cost = local_min_imbalance_cost;
                    best_move = local_min_idx;
                }
                else if (local_min_cut_cost == best_cut_cost && local_min_imbalance_cost < best_imbalance_cost) {
                    best_imbalance_cost = local_min_imbalance_cost;
                    best_move = local_min_idx;
				}
                /*gobal_cost_not_changed_num += local_cost_not_changed_num;
                gobal_out_of_lb_ub_num += local_out_of_lb_ub_num;*/
            }
        }
		cut_distribution.close();
		imbalance_distribution.close();
        //fout2 << "the num of move not change cost: " << gobal_cost_not_changed_num << endl;
        //fout2 << "the num of move out of lb/ub: " << gobal_out_of_lb_ub_num/** 1.0 / num_explores */ << endl;
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
    printf("average time: %.2lf seconds\n", total_time / num_trails);
    return 1;
}
