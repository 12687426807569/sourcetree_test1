using namespace std;
#include<fstream>
#include <iostream>
#include <sstream>

int main()
{
	string file_name_base = "../../ISPD_benchmark/ibm01";
	ifstream fin(file_name_base + ".hgr");
	int cell_num,net_num,fmt;
	string line;
	getline(fin, line);//ignore the first line
	istringstream iss(line);
	iss >>  net_num>>cell_num >>fmt;
	if(!iss.fail()){
		cout << "我還沒做fmt";
		return -1;
	}
	ofstream fout(file_name_base + ".cells");
	for (int i = 1; i <= cell_num; i++)
	{
		fout << 'c' << i << ' ' << 1 << endl;
	}
	fout.close();
	fout.open(file_name_base + ".nets");
	//fin.ignore();//ignore the rest of the line
	int cell_id ;
	while (getline(fin, line))
	{
		istringstream iss(line);
		fout << "NET n1 {";
		while (iss >> cell_id)
			fout << " c" << cell_id;
		fout << " }\n";
	}
	fout.close();
	cout << "finished" << endl;
	return 0;
}