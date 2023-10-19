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
		//保存匹配结果
		double cost = 0;
		for (vector<pair<pair<int, int>, double>>::iterator it = mat.begin();
			it != mat.end(); it++)
			cost += (*it).second;
		Util::setInfo(item, "Greedy", cost, (double)(finish - start) / CLOCKS_PER_SEC, cos.WorkingSetSize / 1024, mat.size());
	}

	vector<pair<pair<int, int>, double>> matching(double budget)
	{
		cout << "matching..." << endl;
		double cost = 0; //已经花费的预算
		vector<pair<pair<int, int>, double>> mat;
		for (int i = 0; i < WORKER_NUM; i++)
		{
			if (Util::greater(cost, budget, 1)) break;
			pair<pair<int, int>, double> p = Util::nearestTask_rt(i + 1, FINF);//在RTree中寻找最近任务
			if (p.first.second != -1 && Util::greater(budget, cost + p.second, 1))//找到最匹配任务
			{
				mat.push_back(p);
				task[p.first.second - 1].matched = 1;
				//从rtree中删除任务节点
				double minpoint[2] = { task[p.first.second - 1].pos_x,task[p.first.second - 1].pos_y };
				tree.Remove(minpoint,minpoint, p.first.second);
				cost = cost + p.second;
				cout << p.first.first << "->" << p.first.second << ":" << p.second <<"total:"<< cost << endl;
			}
		}
		return mat;
	}
};
