#pragma once
#include "Util.h"
class Greedy_RT
{
private:
	double threshold;
	double budget;
public:
	Greedy_RT(double budget,double threshold)
	{
		this->budget = budget;
		this->threshold = threshold;
	}
	Greedy_RT()
	{
	}
	~Greedy_RT()
	{

	}
	//运行一次Greedy-RT
	void run_single(Info &item)
	{
		Util::loadWorkers(wfile.c_str());
		Util::loadTasks(tfile.c_str());
		cout << "building Rtree..." << endl;
		Util::buildRtree();
		PROCESS_MEMORY_COUNTERS cos;
		clock_t start, finish;
		start = clock();
		vector<pair<pair<int, int>, double>> mat = matching();
		finish = clock();
		::GetProcessMemoryInfo(GetCurrentProcess(), &cos, sizeof(cos));

		double cost = 0;
		for (vector<pair<pair<int, int>, double>>::iterator it = mat.begin();
			it != mat.end(); it++)
			cost += (*it).second;
		Util::setInfo(item, "Greedy-RT", cost, (double)(finish - start) / CLOCKS_PER_SEC, cos.WorkingSetSize / 1024, mat.size());
	}
	//运行多次Greedy-RT，求均值
	vector<pair<pair<int, int>, double>> run_ave(double &time,double &memory)
	{
		Util::loadWorkers(wfile.c_str());
		Util::loadTasks(tfile.c_str());
		cout << "building Rtree..." << endl;
		Util::buildRtree();
		PROCESS_MEMORY_COUNTERS cos;
		clock_t start, finish;
		start = clock();
		vector<pair<pair<int, int>, double>> mat = matching();
		finish = clock();
		time+= (double)(finish - start) / CLOCKS_PER_SEC;
		::GetProcessMemoryInfo(GetCurrentProcess(), &cos, sizeof(cos));
		memory += cos.WorkingSetSize / 1024;
		return mat;
	}

	vector<pair<pair<int, int>, double>> matching()
	{
		cout << "matching..." << endl;
		double cost = 0; 
		vector<pair<pair<int, int>, double>> mat;
		for (int i = 0; i < WORKER_NUM; i++)
		{
			if (Util::greater(cost, budget, 1)) break;
			pair<pair<int, int>, double> p = Util::nearestTask_rt(i + 1,threshold);
			if (p.first.second != -1 && Util::greater(budget, cost + p.second, 1))
			{
				mat.push_back(p);
				task[p.first.second - 1].matched = 1;
				//从rtree中删除任务节点
				double minpoint[2] = { task[p.first.second - 1].pos_x,task[p.first.second - 1].pos_y };
				tree.Remove(minpoint, minpoint, p.first.second);
				cost = cost + p.second;
			}
		}
		return mat;
	}
};