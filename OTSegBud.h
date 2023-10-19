#pragma once
#include "Util.h"
class OTSegBud {
private:
	//file
	double budget = 0;
	// adjacency matrix (fill this up)
	int n;  //时间段数
	int time;  //timeline长度
	double threshold;
public:
	OTSegBud(double budget, int n, int time, double threshold)
	{
		this->budget = budget;
		this->n = n;
		this->time = time;
		this->threshold = threshold;
	}
	~OTSegBud()
	{

	}
	//running test case
	void test(Info &item)
	{
		Util::loadWorkers(wfile.c_str());
		Util::loadTasks(tfile.c_str());
		cout << "building Rtree..." << endl;
		Util::buildRtree();
		double used = 0; //记录花费总预算
		double allocated = 0;
		//动态创建数组存放各时间段分配的预算
		double *budPortion = new double[SEGS];
		memset(budPortion, 0, SEGS * sizeof(double));
		if (n != SEGS) { cout << "SegBudOpt中n和SEGS不相同" << endl; return; }
		Util::getOptSegBudget(ofile.c_str(), budPortion, n, time, budget);
		//printf("每个时间段应该分配预算额度\n");
		//for (int i = 0; i < SEGS; i++)
		//{
		//	printf("%f\t", budPortion[i]);
		//	if ((i + 1) % 5 == 0)
		//		printf("\n");
		//}
		PROCESS_MEMORY_COUNTERS cos;
		clock_t start, finish;
		start = clock();
		vector<pair<pair<int, int>, double>> mat = matching(budPortion, SEGS);
		finish = clock();
		::GetProcessMemoryInfo(GetCurrentProcess(), &cos, sizeof(cos));

		double cost = 0;
		for (vector<pair<pair<int, int>, double>>::iterator it = mat.begin();
			it != mat.end(); it++)
			cost += (*it).second;
		Util::setInfo(item, "OTSegBud", cost, (double)(finish - start) / CLOCKS_PER_SEC, cos.WorkingSetSize / 1024, mat.size(),n);
	}

	vector<pair<pair<int, int>, double>> matching(double *budPortion, int n)
	{	
		cout << "matching..." << endl;
		double cost = 0;//记录花费总预算
		double allocated = 0;
		int * allocated_flag = new int[n];
		memset(allocated_flag, 0, n * sizeof(int)); //记录预算是否已分配数组清0
		vector<pair<pair<int, int>, double>> mat;
		for (int i = 0,count=0,mc=0; i < WORKER_NUM; i++)
		{
			//计算工人的所在的时间段k
			int k=worker[i].appear*n/time; //所在的时间段下标
			//count++;
			//将分段预算加入到可用预算中
			if (!allocated_flag[k]) {
				//if (k >= 1)
				//{
				//	printf("第%d个时间段出现了%d工人,匹配%d人，加入预算%f,预算%f,花费%f\n", k, count, mc, budPortion[k-1],allocated, cost);
				//	count = 0;
				//	mc = 0;
				//}
				allocated += budPortion[k];
				allocated_flag[k] = 1;

			}
			if (Util::greater(cost, allocated, 1)) continue;
			pair<pair<int, int>, double> p = Util::nearestTask_rt(i + 1, threshold);
			if (p.first.second != -1 && Util::greater(allocated, cost + p.second, 1))
			{
				mat.push_back(p);
				task[p.first.second - 1].matched = 1;
				//从rtree中删除任务节点
				double minpoint[2] = { task[p.first.second - 1].pos_x,task[p.first.second - 1].pos_y };
				tree.Remove(minpoint, minpoint, p.first.second);
				cost = cost + p.second;
				//mc++;
			}
		}
		delete[] allocated_flag;
		return mat;
	}
};
