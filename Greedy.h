#pragma once
#include "Util.h"
class Greedy {
private:
	double budget;
public:
	Greedy(double budget)
	{
		this->budget = budget;
	}
	~Greedy()
	{

	}

	void run(Info &item)
	{
		Util::loadWorkers(wfile.c_str());
		Util::loadTasks(tfile.c_str());
		cout << "building Rtree..." << endl;
		Util::buildRtree();
		PROCESS_MEMORY_COUNTERS cos;
		clock_t start, finish;
		start = clock();
		vector<pair<pair<int, int>, double>> mat = matching(budget);
		finish = clock();
		::GetProcessMemoryInfo(GetCurrentProcess(), &cos, sizeof(cos));
		//����ƥ����
		double cost = 0;
		for (vector<pair<pair<int, int>, double>>::iterator it = mat.begin();
			it != mat.end(); it++)
			cost += (*it).second;
		Util::setInfo(item, "Greedy", cost, (double)(finish - start) / CLOCKS_PER_SEC, cos.WorkingSetSize / 1024, mat.size());
	}

	vector<pair<pair<int, int>, double>> matching(double budget)
	{
		cout << "matching..." << endl;
		double cost = 0; //�Ѿ����ѵ�Ԥ��
		vector<pair<pair<int, int>, double>> mat;
		for (int i = 0; i < WORKER_NUM; i++)
		{
			if (Util::greater(cost, budget, 1)) break;
			pair<pair<int, int>, double> p = Util::nearestTask_rt(i + 1, FINF);//��RTree��Ѱ���������
			if (p.first.second != -1 && Util::greater(budget, cost + p.second, 1))//�ҵ���ƥ������
			{
				mat.push_back(p);
				task[p.first.second - 1].matched = 1;
				//��rtree��ɾ������ڵ�
				double minpoint[2] = { task[p.first.second - 1].pos_x,task[p.first.second - 1].pos_y };
				tree.Remove(minpoint,minpoint, p.first.second);
				cost = cost + p.second;
				cout << p.first.first << "->" << p.first.second << ":" << p.second <<"total:"<< cost << endl;
			}
		}
		return mat;
	}
};
