#pragma once
#include "Util.h"
class OTSegBud {
private:
	//file
	double budget = 0;
	// adjacency matrix (fill this up)
	int n;  //ʱ�����
	int time;  //timeline����
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
		double used = 0; //��¼������Ԥ��
		double allocated = 0;
		//��̬���������Ÿ�ʱ��η����Ԥ��
		double *budPortion = new double[SEGS];
		memset(budPortion, 0, SEGS * sizeof(double));
		if (n != SEGS) { cout << "SegBudOpt��n��SEGS����ͬ" << endl; return; }
		Util::getOptSegBudget(ofile.c_str(), budPortion, n, time, budget);
		//printf("ÿ��ʱ���Ӧ�÷���Ԥ����\n");
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
		double cost = 0;//��¼������Ԥ��
		double allocated = 0;
		int * allocated_flag = new int[n];
		memset(allocated_flag, 0, n * sizeof(int)); //��¼Ԥ���Ƿ��ѷ���������0
		vector<pair<pair<int, int>, double>> mat;
		for (int i = 0,count=0,mc=0; i < WORKER_NUM; i++)
		{
			//���㹤�˵����ڵ�ʱ���k
			int k=worker[i].appear*n/time; //���ڵ�ʱ����±�
			//count++;
			//���ֶ�Ԥ����뵽����Ԥ����
			if (!allocated_flag[k]) {
				//if (k >= 1)
				//{
				//	printf("��%d��ʱ��γ�����%d����,ƥ��%d�ˣ�����Ԥ��%f,Ԥ��%f,����%f\n", k, count, mc, budPortion[k-1],allocated, cost);
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
				//��rtree��ɾ������ڵ�
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
