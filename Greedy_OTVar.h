#pragma once
#include "Util.h"
class Greedy_OTVar {
private:
	double *threshold;
	double budget;
	int n;
	int time;
public:
	Greedy_OTVar(double budget, double threshold[],int n,int time)
	{
		this->budget = budget;
		this->threshold = threshold;
		this->n = n;
		this->time = time;
	}
	~Greedy_OTVar()
	{

	}
	//test
	void test(Info &item)
	{
		//load data
		Util::loadWorkers(wfile.c_str());
		Util::loadTasks(tfile.c_str());
		cout << "构建任务Rtree中..." << endl;
		Util::buildRtree();
		PROCESS_MEMORY_COUNTERS cos;
		clock_t start, finish;
		start = clock();
		vector<pair<pair<int, int>, double>> mat = matching();
		finish = clock();
		printf("match time:%f\n", (double)(finish - start) / CLOCKS_PER_SEC);
		::GetProcessMemoryInfo(GetCurrentProcess(), &cos, sizeof(cos));
		cout << "memory size:" << cos.WorkingSetSize / 1024 << endl;
		cout << "matching size:" << mat.size() << endl;
		//save statistic information in a .xlsx file
		double cost = 0;
		for (vector<pair<pair<int, int>, double>>::iterator it = mat.begin();
			it != mat.end(); it++)
			cost += (*it).second;
		cout << "cost:" << cost << endl;
		Util::setInfo(item, "Greedy_OTVar", cost, (double)(finish - start) / CLOCKS_PER_SEC, cos.WorkingSetSize / 1024, mat.size(),n);
	}
	//matching 
	vector<pair<pair<int, int>, double>> matching()
	{
		//process the arrive worker one by one
		double cost = 0; //the budget has been used.
		vector<pair<pair<int, int>, double>> mat;
		for (int i = 0; i < WORKER_NUM; i++)
		{
			// worker[i].id, worker[i].pos_x, worker[i].pos_y,worker[i].appear,worker[i].deadline
			int at = worker[i].appear;
			double thr = threshold[at*n / time];
			pair<pair<int, int>, double> p = Util::nearestTask_rt(i + 1, thr); //search the nearest task.
			if (fabs(budget - cost) <= eps) break; //budget has been spend out
			if (p.first.second != -1 && Util::greater(budget, cost + p.second, 1)) //find it
			{
				mat.push_back(p);
				//p.first.second matched, transform logical id to array id
				task[p.first.second - 1].matched = 1;
				//从rtree中删除任务节点
				double minpoint[2] = { task[p.first.second - 1].pos_x,task[p.first.second - 1].pos_y };
				tree.Remove(minpoint, minpoint, p.first.second);
				cost = cost + p.second;
				//cout <<"出现时刻："<<at<<" 时间段："<<at*n/time<<" 阈值："<<thr<< " w" << p.first.first << "->" << "t" << p.first.second << "cost:" << p.second << endl;
			}
		}
		return mat;
	}

};


