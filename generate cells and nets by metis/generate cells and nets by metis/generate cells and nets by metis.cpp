using namespace std;
#include<fstream>
#include <iostream>
#include <sstream>

int main()
{
	string file_name_base = "../../simple_testcase/simple_graph_testcase";
	ifstream fin(file_name_base + ".metis");
	ofstream fout(file_name_base + ".cells");
	int cell_num;
	fin >> cell_num;;
	for (int i = 1; i <= cell_num; i++)
	{
		fout << 'c' << i << ' ' << 1 << endl;
	}
	fout.close();
	fout.open(file_name_base + ".nets");
	int net_num;
	fin >> net_num;
	fin.ignore();//ignore the rest of the line
	int cell_id = 1;
	string line;
	while (getline(fin, line))
	{
		istringstream iss(line);
		int neighbor;
		while (iss >> neighbor)
			if (neighbor > cell_id) //不接受自己連自己，也不接受重複邊(無向圖(1,2)=(2,1))
				fout << "NET n1 { c" << cell_id << " c" << neighbor << " }\n";
		cell_id++;
	}
	fout.close();
	cout << "finished" << endl;
	return 0;
}