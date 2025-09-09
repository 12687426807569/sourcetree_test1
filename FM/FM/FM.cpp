using namespace std;
#include <iostream>
#include<fstream>
#include<tuple>
#include <chrono>
#include<unordered_map>
#include<unordered_set>
#include <vector>
#include <string>
#include<climits>
#include<set>
#include<array>
#include <list>
#include<cmath>
#define locked true
#define unlocked false
using cell_program_id_type = int;
using cell_input_id_type = int;
using net_id_type = int;
constexpr int A = 0, B = 1;
using net_information = array< unordered_set<cell_program_id_type>, 2>;//if all cells in A.cells_classified_by_partition[B] still exist
const array<char, 2> partition_input_name_s = { 'A','B' };
template <typename Key, typename Value>
struct myMap {
public:
	std::unordered_map<Key, Value> umap;  // O(1) 存取
	std::set<Key> keyset;                 // O(log n) 保序

	// 包裝迭代器
	struct iterator {
	public:
		typename std::set<Key>::iterator sit;
		std::unordered_map<Key, Value>* ump;

		iterator(typename std::set<Key>::iterator s, std::unordered_map<Key, Value>* u) : sit(s), ump(u) {}

		iterator& operator++() { ++sit; return *this; }
		iterator operator++(int) { auto tmp = *this; ++(*this); return tmp; }
		iterator& operator--() { --sit; return *this; }
		bool operator==(const iterator& other) const { return sit == other.sit; }
		bool operator!=(const iterator& other) const { return sit != other.sit; }

		// 直接回傳 std::pair
		std::pair<const Key, Value&> operator*() const {
			return { *sit, (*ump)[*sit] };
		}
		typename unordered_map<Key, Value>::iterator operator->() {//a->b=(a->())->
			typename unordered_map<Key, Value>::iterator it = ump->find(*sit);
			return it;
		}
	};

	// 插入或修改
	// 刪除
	void erase(const Key k) {
		keyset.erase(k);
		umap.erase(k);
	}
	void erase( iterator it) {
		auto k = *(it.sit);
		keyset.erase(it.sit);
		umap.erase(k);
	}

	// operator[]
	typename Value& operator[](const Key& k) {
		if (umap.find(k) == umap.end())
			keyset.insert(k);
		return umap[k];
	}

	// 遍歷 (用自定義 iterator)
	iterator begin() { return iterator(keyset.begin(), &umap); }
	iterator end() { return iterator(keyset.end(), &umap); }

	size_t size() const { return umap.size(); }
	bool empty() const { return umap.empty(); }
};
struct cell_information
{
	int in_partition_name = A;
};
struct partition_information
{
	unordered_set<cell_program_id_type> cells;
	int part_area = 0;
};
struct myGraph {
	vector<int>cell_program_id_to_area;//V[cell_program_id]=area
	vector<cell_input_id_type>cell_program_id_to_input_id;//v[cell_program_id]=cell_input_id
	unordered_map<cell_input_id_type, cell_program_id_type>cell_input_id_to_program_id;//v[cell_input_id]=cell_program_id
	vector<vector<cell_program_id_type>> net_id_to_cell_program_id_s;//v[net_id]={cell_program_id:cell_program_id conneted net_id}
	vector< vector<net_id_type>>cell_program_id_to_net_id_s;//v[cell_program_id]={net_id:net_id conneted cell_program_id}
};
myGraph file_read(const string cell_file, const string net_file)
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
		G.cell_program_id_to_net_id_s.push_back({});
		i++;
	}
	fin.close();
	fin.open(net_file);
	string t,net_name;
	while (fin >> t)
	{
		fin >> net_name >> t;
		vector<cell_program_id_type> cell_id_s;
		while (fin >> t)
		{
			if (t == "}")break;
			cell_program_id_type cell_program_id = G.cell_input_id_to_program_id[stoi(t.substr(1))];
			cell_id_s.push_back(cell_program_id);
			G.cell_program_id_to_net_id_s[cell_program_id].push_back(G.net_id_to_cell_program_id_s.size());
		}
		G.net_id_to_cell_program_id_s.push_back(move(cell_id_s));
	}
	fin.close();
	return G;
}
pair<int, int> my_BALANCE_CRITERION(const vector<int>& cell_program_id_to_area, const double r, const double epsilon) {
	int sum = 0, max_area = 0;
	for (int area : cell_program_id_to_area) {// C++17 structured binding
		sum += area;
		if (area > max_area) max_area = area;
	}
	int  S_opt = ceil(r * sum);//definition of https://chriswalshaw.co.uk/partition/
	int max_S_p =floor( (1 + epsilon) * S_opt);//0.52*sum_area
	if(ceil(r * sum) + max_area> max_S_p) cout<<"imbalance rate too small",exit(1);//min(0.48*sum_area,0.5*sum_area-max_area)
	return {sum - max_S_p,  max_S_p };
}
array< partition_information, 2> PARTITION(const int mean, const myGraph& G) {
	array<partition_information, 2> partition_informations;
	int i = 0, size = G.cell_program_id_to_area.size();
	for (; i < size; ++i) {
		if (partition_informations[A].part_area >= mean)break;
		partition_informations[A].part_area += G.cell_program_id_to_area[i];
		partition_informations[A].cells.insert(i);
	}
	for (; i < size; ++i) {
		partition_informations[B].part_area += G.cell_program_id_to_area[i];
		partition_informations[B].cells.insert(i);
	}
	return partition_informations;
}
pair< vector<net_information>, vector< cell_information> > using_PARTITION(
	const array< partition_information, 2>& partition_informations, const
	myGraph& G) {
	vector<net_information> net_info_temporary;
	vector<  cell_information> cell_info(G.cell_program_id_to_area.size());
	for (int cell_id : partition_informations[A].cells)
		cell_info[cell_id].in_partition_name = A;
	for (int cell_id : partition_informations.at(B).cells)
		cell_info[cell_id].in_partition_name = B;
	for (auto& net : G.net_id_to_cell_program_id_s) {
		net_info_temporary.push_back(net_information());
		auto& a = net_info_temporary.back();
		for (int cell_id : net)
			a[cell_info[cell_id].in_partition_name].insert(cell_id);
	}
	//return { net_info_temporary, cell_info };//call copy constructor,not claim return value in caller frame
	return { move(net_info_temporary), move(cell_info) };//call move constructor
}
vector<int> FS_minus_TE(const myGraph& G, const vector<net_information>& net_info_s) {
	//genTestcase.rb maybe generate a  cell that no net connects it,so we need to check
	vector< int> delta_g_i(G.cell_program_id_to_area.size(), 0);
	for (int partition_name = 0; partition_name < 2; partition_name++) {
		const int another_part_name = A + B - partition_name;
		for (const net_information& net_info_1 : net_info_s) {
			const auto& cell_id_set = net_info_1[partition_name];
			const auto& cell_id_set_another = net_info_1[another_part_name];
			if (cell_id_set.size() == 1)
				delta_g_i[*cell_id_set.begin()] += 1;
			if (cell_id_set.size() == 0)
				for (int cell : cell_id_set_another)
					delta_g_i[cell] -= 1;
		}
	}
	return delta_g_i;//NRVO,不用 move
}
int MAX_GAIN(myMap<int, unordered_set<cell_program_id_type>>& free_cells_sorted_by_delta_g, const int lb,
	const int ub, const array< partition_information, 2>&
	partition_informations_s, const vector<cell_information>&
	cell_info_s, const myGraph& G) {
	for (auto it = free_cells_sorted_by_delta_g.end();
		it != free_cells_sorted_by_delta_g.begin(); ) {
		--it;
		vector<pair<cell_program_id_type, int>> candidate;//cell_id,partition_area*2 after move
		for (auto it1 = it->second.end(); it1 != it->second.begin(); ) {
			it1--;
			int from = cell_info_s.at(*it1).in_partition_name;
			if (partition_informations_s.at(from).part_area - G.cell_program_id_to_area.
				at(*it1) >= lb) {
				candidate.emplace_back(*it1, (partition_informations_s.at(from)
					.part_area - G.cell_program_id_to_area.at(*it1))*2); // Found a valid cell to move	
			}
		}
		if (candidate.empty()) continue;
		const int mid_2 = (ub + lb);
		auto min_cell_id = candidate[0].first;
		int min_num = labs( candidate[0].second - mid_2);
		for (auto& [candidate_it, partition_area_2] : candidate) {
			if (labs(partition_area_2 - mid_2) < min_num) {
				min_num = labs(partition_area_2 - mid_2);
				min_cell_id = candidate_it;
			}
		}
		it->second.erase(min_cell_id);
		if (it->second.empty())
			free_cells_sorted_by_delta_g.erase(it);
		return min_cell_id; // Found a valid cell to move
	}
	return -1;
}
vector<net_id_type> CRITIAL_NETS(const cell_program_id_type cell_id, const myGraph& G, const
	vector<net_information>& net_info_s, const vector<cell_information>
	& cell_info_s) {
	vector<net_id_type>ret;
	int from = cell_info_s.at(cell_id).in_partition_name, to = A + B - from;
	//genTestcase.rb maybe generate a  cell that no net connects it,so we need to check
	for (net_id_type net : G.cell_program_id_to_net_id_s.at(cell_id)) {
		auto& a = net_info_s[net];
		if (a.at(from).size() <= 2 || a.at(to).size() <= 1)
			ret.push_back(net);
	}
	return ret;//NRVO,不用 move
}
void TRY_MOVE(array< partition_information, 2>& partition_info, vector<
	cell_information>& cell_info, vector<net_information>& net_info,
	const cell_program_id_type cell, const char from, const char to, const myGraph& G) {
	cell_info[cell].in_partition_name = to;
	partition_info[from].part_area -= G.cell_program_id_to_area.at(cell);
	partition_info[from].cells.erase(cell);
	partition_info[to].part_area += G.cell_program_id_to_area.at(cell);
	partition_info[to].cells.insert(cell);
	for (int net : G.cell_program_id_to_net_id_s.at(cell)) {
		auto& a = net_info[net];
		a[from].erase(cell);
		a[to].insert(cell);
	}
}
void UPDATE_GAIN(vector< int>& delta_g_i, myMap<int, unordered_set<int>>&
	free_cells_sorted_by_delta_g, vector< cell_information>& cell_info_temporary,
	const char from, const char to, const net_information& cells_classified_by_partition, const
	vector< char>& status) {
	if (cells_classified_by_partition.at(to).size() == 2)
		for (int cell_id : cells_classified_by_partition.at(to))
			if (status.at(cell_id) == unlocked) {
				free_cells_sorted_by_delta_g[delta_g_i[cell_id]].erase(
					cell_id);
				if (free_cells_sorted_by_delta_g[delta_g_i[cell_id]].empty())
					free_cells_sorted_by_delta_g.erase(delta_g_i[cell_id]);
				delta_g_i[cell_id]--;
				free_cells_sorted_by_delta_g[delta_g_i[cell_id]].insert(cell_id);
			}
	if (cells_classified_by_partition.at(to).size() == 1)
		for (int cell_id : cells_classified_by_partition.at(from))
			if (status.at(cell_id) == unlocked) {
				free_cells_sorted_by_delta_g[delta_g_i[cell_id]].erase(cell_id);
				if (free_cells_sorted_by_delta_g[delta_g_i[cell_id]].empty())
					free_cells_sorted_by_delta_g.erase(delta_g_i[cell_id]);
				delta_g_i[cell_id]++;
				free_cells_sorted_by_delta_g[delta_g_i[cell_id]].insert( cell_id);
			}
	if (cells_classified_by_partition.at(from).size() == 0)
		for (int cell_id : cells_classified_by_partition.at(to))
			if (status.at(cell_id) == unlocked) {
				free_cells_sorted_by_delta_g[delta_g_i[cell_id]].erase(cell_id);
				if (free_cells_sorted_by_delta_g[delta_g_i[cell_id]].empty())
					free_cells_sorted_by_delta_g.erase(delta_g_i[cell_id]);
				delta_g_i[cell_id]--;
				free_cells_sorted_by_delta_g[delta_g_i[cell_id]].insert( cell_id);
			}
	if (cells_classified_by_partition.at(from).size() == 1)
		for (int cell_id : cells_classified_by_partition.at(from))
			if (status.at(cell_id) == unlocked) {
				free_cells_sorted_by_delta_g[delta_g_i[cell_id]].erase(cell_id);
				if (free_cells_sorted_by_delta_g[delta_g_i[cell_id]].empty())
					free_cells_sorted_by_delta_g.erase(delta_g_i[cell_id]);
				delta_g_i[cell_id]++;
				free_cells_sorted_by_delta_g[delta_g_i[cell_id]].insert(cell_id);
			}
}
pair<int, int> BEST_MOVES(const vector<tuple<int, cell_program_id_type, int>>&
	order) {
	int Gm = INT_MIN, m = 0, sum = 0;
	for (int i = 0; i < order.size();) {
		sum += get<2>(order[i]);
		i++;
		if (Gm < sum) {
			Gm = sum;
			m = i;
		}
	}
	return { Gm, m };
}
void CONFIRM_MOVES(array< partition_information, 2>&
	partition_informations, const vector<tuple<int, cell_program_id_type, int>>&
	order, const int m, const myGraph& G) {/*
	for (int i = 0; i < m; i++) {
		char from = get<0>(order[i]), to = A + B - from;
		cell_program_id_type cell_id = get<1>(order[i]);
		partition_informations[from].part_area -= G.cell_program_id_to_area.at(cell_id);
		partition_informations[from].cells.erase(cell_id);
		partition_informations[to].part_area += G.cell_program_id_to_area.at(cell_id);
		partition_informations[to].cells.insert(cell_id);
	}*/
	for(int i=m;i<order.size();i++) {
		auto from = get<0>(order[i]), to = A + B - from;
		cell_program_id_type cell_id = get<1>(order[i]);
		partition_informations[from].part_area += G.cell_program_id_to_area.at(cell_id);
		partition_informations[from].cells.insert(cell_id);
		partition_informations[to].part_area -= G.cell_program_id_to_area.at(cell_id);
		partition_informations[to].cells.erase(cell_id);
	}
}
#define TEST_CASE_NUM 1
int main()
{
	for (int test_case = 1; test_case <= TEST_CASE_NUM; test_case++) {
		chrono::high_resolution_clock::time_point start = chrono::high_resolution_clock::now();
		string file_name_base = /*string("../../colab_testcase 10000 100 10000 2/test") + to_string(test_case);
			//*/"../../testcase/delaunay_n10";
		string cell_file = file_name_base + ".cells";
		string net_file = file_name_base + ".nets";
		ofstream fout(file_name_base + ".partitions");
		const myGraph G = file_read(cell_file, net_file);
		int total_area = 0;
		for (int area : G.cell_program_id_to_area)
			total_area += area;
		const int ub = floor(0.5 * total_area + total_area / 20.0),lb=total_area - ub ;//definition of github.com/EricLu1218/Physical_Design_Automation/blob/main/Two-way_Min-cut_Partitioning/CS613500_HW2_spec.pdf		
		//const auto [lb, ub] = my_BALANCE_CRITERION(G.cell_program_id_to_area, 0.5, 0.04);	//c++17 auto[]
		array< partition_information, 2> partition_informations = PARTITION((ub + lb) / 2, G);
		int Gm = INT_MAX;
		vector<tuple<int, cell_program_id_type, int>> order;//(from,cell,delta_g)
		while (Gm > 0) {
			order.clear();//capacity不變
			auto [net_info_temporary, cell_info_temporary] = using_PARTITION(
				partition_informations, G);
			int cost = 0;
			for (auto& a : net_info_temporary) {
				if (a[A].size() > 0 && a[B].size() > 0) {
					cost++;
				}
			}
			cout << "cost: " << cost << endl;
			//改try_move變真的movr,comfirm move把不需要的move取消，省暫時partition的複製時間
			//array< partition_information, 2> partition_informations_temporary =partition_informations;
			vector< int> delta_g_i = FS_minus_TE(G, net_info_temporary);
			myMap<int, unordered_set<cell_program_id_type>> free_cells_sorted_by_delta_g;//可用ordered map+存key的set優化 map 不新增key的[]操作
			vector<char/*bool*/>status(G.cell_program_id_to_area.size());//vector<bool>存取慢
			for (int cell_id = 0, size = G.cell_program_id_to_area.size(); cell_id < size; cell_id++) {
				free_cells_sorted_by_delta_g[delta_g_i[cell_id]].insert(cell_id);
				status[cell_id] = unlocked;
			}
			while (free_cells_sorted_by_delta_g.size() != 0) {
				int cell = MAX_GAIN(free_cells_sorted_by_delta_g, lb, ub,
					partition_informations, cell_info_temporary, G);
				order.emplace_back(cell_info_temporary[cell].in_partition_name,
					cell, delta_g_i[cell]);
				vector<net_id_type>critical_nets = CRITIAL_NETS(cell, G, net_info_temporary,
					cell_info_temporary);
				auto from = cell_info_temporary[cell].in_partition_name, to = A + B - from;
				TRY_MOVE(partition_informations, cell_info_temporary,
					net_info_temporary, cell, from, to, G);//move
				status[cell] = locked;
				for (int net : critical_nets)
					UPDATE_GAIN(delta_g_i, free_cells_sorted_by_delta_g, cell_info_temporary, from, to,
						net_info_temporary[net], status);
			}
			int m; tie(Gm, m) = BEST_MOVES(order);
			if (Gm > 0)
				CONFIRM_MOVES(partition_informations, order, m, G);//order[0] to order[m-1] are the best moves,cancel m to order.size()-1 moves
		}
		auto [net_info_temporary, cell_info_temporary] = using_PARTITION(
			partition_informations, G);
		fout << "cut_size ";
		int cost = 0;
		for (auto& a : net_info_temporary)
			if (a[A].size() > 0 && a[B].size() > 0)
				cost++;
		cout << "test case " << test_case << ",cost= " << cost << endl;
		fout << cost << endl;
		for (int i = 0; i < 2; i++) {
			auto& partition_info = partition_informations[i];
			fout << partition_input_name_s[i] << " " << partition_info.cells.size() << endl;
			for (int cell_id : partition_info.cells) {
				fout << "C" << G.cell_program_id_to_input_id[cell_id] << endl;
			}
		}
		fout.close();
		cout << "Time taken: "
			<< chrono::duration_cast<chrono::milliseconds>(chrono::high_resolution_clock::now() - start).count()
			<< " milliseconds" << endl;
	}
	cout << "finished" << endl;
	return 0;
}