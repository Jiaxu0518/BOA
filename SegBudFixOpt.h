#pragma once
/*
set:NN
call:
Opt *opt=new Opt(10);    create by new method is demanded
opt->test();
delete opt;
*
*/
#include "Util.h"
#include "psapi.h"
#define CLR(a, x) memset( a, x, sizeof( a ) )
#define BUBL { \
    t = q[i]; q[i] = q[j]; q[j] = t; \
    t = inq[q[i]]; inq[q[i]] = inq[q[j]]; inq[q[j]] = t; }

// Dijkstra's using non-negative edge weights (cost + potential)
#define Pot(u,v) (d[u] + pi[u] - pi[v])
#pragma comment(lib, "Psapi.lib")
using namespace std;
// the maximum number of vertices + 1

class SegBudFixOpt {
private:
	//file
	double budget = 0;
	// adjacency matrix (fill this up)
	uint8_t cap[NN][NN];
	//float cost[NN][NN];
	//cost matrix
	float s[WORKER_NUM];
	float e[TASK_NUM];
	float cross[WORKER_NUM][TASK_NUM];
	// flow network and adjacency list
	// cost per unit of flow matrix (fill this up)
	uint8_t fnet[NN][NN];
	short int adj[NN][NN], deg[NN];
	// Dijkstra's predecessor, depth and priority queue
	short int par[NN], q[NN], inq[NN], qs;
	float d[NN];
	// Labelling function
	float pi[NN];
	//统计各时间片可用的任务数和工人数
	vector<int> avail_workers;
	vector<int> avail_tasks;
	int n;  //时间段数
	int time;  //timeline长度
	double threshold;
public:
	SegBudFixOpt(double budget, int n, int time, double threshold)
	{
		this->budget = budget;
		this->n = n;
		this->time = time;
		this->threshold = threshold;
	}
	~SegBudFixOpt()
	{

	}
	float cost(int i, int j) {
		if (i >= 1 && i <= WORKER_NUM &&j >= WORKER_NUM + 1 && j <= WORKER_NUM + TASK_NUM) return cross[i - 1][j - WORKER_NUM - 1];
		else if (j >= 1 && j <= WORKER_NUM &&i >= WORKER_NUM + 1 && i <= WORKER_NUM + TASK_NUM) return cross[i - WORKER_NUM - 1][j - 1];
		else if (j >= 1 && j <= WORKER_NUM &&i == 0) return s[j - 1];
		else if (i >= 1 && i <= WORKER_NUM &&j == 0) return s[i - 1];
		else if (i == WORKER_NUM + TASK_NUM + 1 && j >= WORKER_NUM + 1 && j <= WORKER_NUM + TASK_NUM) return e[j - WORKER_NUM - 1];
		else if (j == WORKER_NUM + TASK_NUM + 1 && i >= WORKER_NUM + 1 && i <= WORKER_NUM + TASK_NUM) return e[i - WORKER_NUM - 1];
		else return FINF;
	}
	//running test case
	void test(Info &item)
	{
		//匹配方案
		vector<pair<pair<int, int>, double>> mat;
		//load workers' and tasks' information from .xlsx files.
		Util::loadWorkers(wfile.c_str());
		Util::loadTasks(tfile.c_str());
		double used = 0; //记录花费总预算
		double allocated = 0;
		//动态创建数组存放各时间段分配的预算
		double *budPortion = new double[SEGS];
		memset(budPortion, 0, SEGS * sizeof(double));
		if (n != SEGS) { cout << "SegBudOpt中n和SEGS不相同" << endl; return; }
		Util::getOptSegBudget(ofile.c_str(), budPortion, n, time, budget);
		PROCESS_MEMORY_COUNTERS cos;
		clock_t start, finish;
		start = clock();
		cout << "matching... total segs " << n;
		for (int i = 1; i <= n; i++) //每一段统计一次
		{
			//cout << "matching for the " << i << " segment...." << endl;
			float fcost = 0;
			allocated += budPortion[i - 1]; //每个时间段加入该段的预算
											//可以用工人和任务
			searchSEGAvail_r(i);
			//cout << "第" << i << "时间段可用工人数：" << avail_workers.size() << "可用任务数：" << avail_tasks.size() << endl;
			//清空数组
			clearArrays();
			//匹配任务
			// numV is the number of vertices,including source and sink vertex.
			int numV = WORKER_NUM + TASK_NUM + 2;
			memset(cap, 0, sizeof(cap));//clear 0
										//specify the id of source and sink vertex 
			int s = 0, t = numV - 1;
			//cout << "nWorkers:" << WORKER_NUM << " nTask:" << TASK_NUM << " nVertices:" << numV << " sid:" << s << " tid:" << t << endl;
			//build bipartition graph
			buildGraph(WORKER_NUM, TASK_NUM);
			//matching
			//cout << "mathcing..." << endl;
			int flow = mcmf4(numV, s, t, fcost, allocated - used); //allocated - used是该轮分配的预算
			used += fcost;//把该段消耗掉的预算计入到总花费中
						  //提取最优匹配并加入结果集
			vector<pair<pair<int, int>, double>> m = optMatching();
			if (m.size() != 0)
				mat.insert(mat.end(), m.begin(), m.end());
			//删除第i个时间段执行过后，已经过期和匹配的工人任务
			removeSEGAvail(m);
			//cout << "第" << i << "时间段新增匹配数:" << flow
			//	<< " 该轮投入的预算：" << allocated - used + fcost
			//	<< " 用掉的预算：" << fcost
			//	<< " 总花费：" << used
			//	<< " 总的匹配数量：" << mat.size() << endl;
		}
		finish = clock();
		printf("match time:%f\n", (double)(finish - start) / CLOCKS_PER_SEC);
		::GetProcessMemoryInfo(GetCurrentProcess(), &cos, sizeof(cos));
		//save statistic information in a .xlsx file
		Util::setInfo(item, "SegBudFixOpt", used, (double)(finish - start) / CLOCKS_PER_SEC, cos.WorkingSetSize / 1024, mat.size(), n);
		cout << "memory size:" << cos.WorkingSetSize / 1024 << endl;
		delete[] budPortion;
	}
	bool dijkstra(int n, int s, int t)
	{
		for (int i = 0; i < NN; i++)
			d[i] = 0x3f;
		//CLR(d, 0x3F);
		CLR(par, -1);
		CLR(inq, -1);
		d[s] = qs = 0;
		inq[q[qs++] = s] = 0;
		par[s] = n;
		while (qs)
		{
			// get the minimum from q and bubble down
			int u = q[0]; inq[u] = -1;
			q[0] = q[--qs];
			if (qs) inq[q[0]] = 0;
			for (int i = 0, j = 2 * i + 1, t; j < qs; i = j, j = 2 * i + 1)
			{
				if (j + 1 < qs && Util::greater(d[q[j]], d[q[j + 1]], 0)) j++;
				if (Util::greater(d[q[j]], d[q[i]], 1)) break;
				BUBL;
			}
			// relax edge (u,i) or (i,u) for all i;
			for (int k = 0, v = adj[u][k]; k < deg[u]; v = adj[u][++k])
			{
				// try undoing edge v->u
				if (fnet[v][u] && Util::greater(d[v], Pot(u, v) - cost(v, u), 0))
					d[v] = Pot(u, v) - cost(v, par[v] = u);

				// try using edge u->v
				if (fnet[u][v] < cap[u][v] && Util::greater(d[v], Pot(u, v) + cost(u, v), 0))
					d[v] = Pot(u, v) + cost(par[v] = u, v);

				if (par[v] == u)
				{
					// bubble up or decrease key
					if (inq[v] < 0) { inq[q[qs] = v] = qs; qs++; }
					for (int i = inq[v], j = (i - 1) / 2, t;
						Util::greater(d[q[j]], d[q[i]], 0); i = j, j = (i - 1) / 2)
						BUBL;
				}
			}
		}
		for (int i = 0; i < n; i++) if (Util::greater(FINF, pi[i], 0)) pi[i] += d[i];
		return par[t] >= 0;
	}
	//min-cost-max-flow
	int mcmf4(int n, int s, int t, float &fcost, double budget)
	{
		// build the adjacency list
		CLR(deg, 0);
		for (int i = 0; i < n; i++)
			for (int j = 0; j < n; j++)
				if (cap[i][j] || cap[j][i]) adj[i][deg[i]++] = j;

		CLR(fnet, 0);
		CLR(pi, 0);
		int flow = 0;
		// repeatedly, find a cheapest path from s to t
		while (dijkstra(n, s, t))
		{
			// get the bottleneck capacity
			int bot = INT_MAX;
			for (int v = t, u = par[v]; v != s; u = par[v = u])
			{
				int tem = fnet[v][u] ? fnet[v][u] : (cap[u][v] - fnet[u][v]);
				bot = bot<tem ? bot : tem;
			}
			// if fcost beyonds budget when adding the current flow, give up and terminate
			float nextb = fcost;
			for (int v = t, u = par[v]; v != s; u = par[v = u])//routine of the last flow
				if (fnet[v][u])
				{
					nextb -= bot * cost(v, u);
				}
				else
				{
					nextb += bot * cost(u, v);
				}
			if (Util::greater(nextb, budget, 0))
				break; //teminate in advance
					   // update the flow network
			for (int v = t, u = par[v]; v != s; u = par[v = u])
				if (fnet[v][u])
				{
					fnet[v][u] -= bot; fcost -= bot * cost(v, u);
				}
				else
				{
					fnet[u][v] += bot; fcost += bot * cost(u, v);
				}
			flow += bot;
		}
		return flow;
	}

	//build graph with workers and tasks nodes
	void buildGraph(int workerNum, int taskNum)
	{
		//start index 0, terminal index workerNum + taskNum + 1, number of vertexs is (|workers|+|tasks|)
		// add edge between s and worker vertex
		//cout << "构建空图（容量设置为0，travel cost设置为FINF）" << endl;
		//构建空图（容量设置为0，travel cost设置为FINF）
		for (int i = 1; i <= workerNum; i++) {
			s[i - 1] = FINF; // cost[b][a] = c;
			cap[0][i] = 0; // cap[b][a] = cp;
		}
		// add edge between t and task vertex
		for (int i = workerNum + 1; i <= workerNum + taskNum; i++) {
			e[i - workerNum - 1] = FINF; // cost[b][a] = c;
			cap[i][workerNum + taskNum + 1] = 0; // cap[b][a] = cp;
		}
		// add edge between worker and task vertices
		for (int taskid = workerNum + 1; taskid <= workerNum + taskNum; taskid++) {
			for (int workerid = 1; workerid <= workerNum; workerid++)
			{
				//transform the id in graph to the logical id
				int wid = workerid;
				int tid = taskid - workerNum;
				cross[workerid - 1][taskid - workerNum - 1] = FINF; // cost[b][a] = c;
				cap[workerid][taskid] = 0; // cap[b][a] = cp;
			}
		}
		//cout << "在空图基础上构建新图" << endl;
		//在空图基础上构建新图
		for (vector<int>::iterator itw = avail_workers.begin(); itw != avail_workers.end(); itw++)
		{
			int i = (*itw);//获得工人的逻辑id
			s[i - 1] = 0; // cost[b][a] = c;
			cap[0][i] = 1; // cap[b][a] = cp;
		}
		// add edge between t and task vertex
		for (vector<int>::iterator itt = avail_tasks.begin(); itt != avail_tasks.end(); itt++)
		{
			int i = workerNum + (*itt);
			e[i - workerNum - 1] = 0; // cost[b][a] = c;
			cap[i][workerNum + taskNum + 1] = 1; // cap[b][a] = cp;
		}
		// add edge between worker and task vertices
		for (vector<int>::iterator itt = avail_tasks.begin(); itt != avail_tasks.end(); itt++)
			for (vector<int>::iterator itw = avail_workers.begin(); itw != avail_workers.end(); itw++)
			{
				int wid = (*itw); //逻辑id
				int tid = (*itt);
				int workerid = wid; //在图中节点id
				int taskid = tid + workerNum;
				double dis;
				if (mode)
					dis = Util::Manhattan(wid, tid);  //syn
				else
					dis = Util::Euclidean_Real(wid, tid); //real
				if (fabs(dis + 1)>eps && Util::greater(threshold, dis, 1)) //worker can arrive before deadline
				{
					cross[workerid - 1][taskid - workerNum - 1] = dis; // cost[b][a] = c;
					cap[workerid][taskid] = 1; // cap[b][a] = cp;
				}
				else
				{
					cross[workerid - 1][taskid - workerNum - 1] = FINF;
					cap[workerid][taskid] = 0; // cap[b][a] = cp;
				}
			}
	}
	//return optimal matching pairs from fnet[NN][NN],where worker in[1,WORKER_NUM],task in[WORKER_NUM+1,WORKER_NUM+TASK_NUM]
	vector<pair<pair<int, int>, double>> optMatching()
	{
		vector<pair<pair<int, int>, double>> v;
		for (int windex = 1; windex <= WORKER_NUM; windex++)
			for (int tindex = WORKER_NUM + 1; tindex <= WORKER_NUM + TASK_NUM; tindex++)
			{
				if (fnet[windex][tindex] == 1)
					v.push_back(make_pair(make_pair(windex, tindex - WORKER_NUM), cost(windex, tindex))); //保存的是工人和任务的逻辑id
			}
		return v;
	}
	//寻找某时间段可用的工人和任务（每个工人和任务随时退出）
	void searchSEGAvail_r(int sid) //逻辑时间段id
	{
		//增加可用工人
		for (int i = 0; i < WORKER_NUM; i++) {
			if (Util::greater(worker[i].appear, 1.0*(sid - 1)*time / n, 1) && Util::greater(1.0*sid*time / n, worker[i].appear, 0))
			{
				avail_workers.push_back(i + 1);//加入工人的逻辑ID
			}
		}
		//增加可用任务
		for (int i = 0; i < TASK_NUM; i++) {
			if (Util::greater(task[i].release, 1.0*(sid - 1)*time / n, 1) && Util::greater(1.0*sid*time / n, task[i].release, 0))
			{
				avail_tasks.push_back(i + 1);//加入任务的逻辑ID
			}
		}
		//删除已经过期的工人,在第sid个时间段截止时间过期的工人
		for (vector<int>::iterator it = avail_workers.begin(); it != avail_workers.end(); )
		{
			int wid = *it;
			if (Util::greater(1.0*(sid-1)*time / n, worker[wid - 1].deadline, 0)) //工人过期则删除
				avail_workers.erase(it);
			else
				it++;
		}
		//删除已经过期的任务，在sid个时间段不会出现的工人
		for (vector<int>::iterator it = avail_tasks.begin(); it != avail_tasks.end(); )
		{
			int tid = *it;
			if (Util::greater(1.0*(sid-1)*time / n, task[tid - 1].deadline, 0)) //任务过期则删除
				avail_tasks.erase(it);
			else
				it++;
		}
	}
	//寻找某时间段可用的工人和任务（每个工人和任务即使deadline到了也需要等到匹配完成后才允许退出）
	void searchSEGAvail(int sid) //逻辑时间段id
	{
		//删除已经过期的工人,在第sid个时间段不会出现的工人
		for (vector<int>::iterator it = avail_workers.begin(); it != avail_workers.end(); )
		{
			int wid = *it;
			if (Util::greater(1.0*(sid - 1)*time / n, worker[wid - 1].deadline, 0)) //工人过期则删除
				avail_workers.erase(it);
			else
				it++;
		}
		//删除已经过期的任务，在sid个时间段不会出现的工人
		for (vector<int>::iterator it = avail_tasks.begin(); it != avail_tasks.end(); )
		{
			int tid = *it;
			if (Util::greater(1.0*(sid - 1)*time / n, task[tid - 1].deadline, 0)) //任务过期则删除
				avail_tasks.erase(it);
			else
				it++;
		}
		//增加可用工人
		for (int i = 0; i < WORKER_NUM; i++) {
			if (Util::greater(worker[i].appear, 1.0*(sid - 1)*time / n, 1) && Util::greater(1.0*sid*time / n, worker[i].appear, 0))
			{
				avail_workers.push_back(i + 1);//加入工人的逻辑ID
			}
		}
		//增加可用任务
		for (int i = 0; i < TASK_NUM; i++) {
			if (Util::greater(task[i].release, 1.0*(sid - 1)*time / n, 1) && Util::greater(1.0*sid*time / n, task[i].release, 0))
			{
				avail_tasks.push_back(i + 1);//加入任务的逻辑ID
			}
		}
	}
	//删除某时间段执行过后，已经匹配的工人任务
	void removeSEGAvail(vector<pair<pair<int, int>, double>> &m)
	{
		//删除已经匹配的工人和任务
		for (vector<pair<pair<int, int>, double>>::iterator it = m.begin(); it != m.end(); it++)
		{
			pair<pair<int, int>, double> p = *it;
			vector<int>::iterator itw = find(avail_workers.begin(), avail_workers.end(), p.first.first); //查找并删除工人
			if (itw != avail_workers.end())
				avail_workers.erase(itw);
			vector<int>::iterator itt = find(avail_tasks.begin(), avail_tasks.end(), p.first.second); //查找并删除任务
			if (itt != avail_tasks.end())
				avail_tasks.erase(itt);
		}
	}
	void clearArrays()
	{
		memset(cap, 0, NN*NN * sizeof(uint8_t));
		memset(s, 0, WORKER_NUM * sizeof(float));
		memset(e, 0, TASK_NUM * sizeof(float));
		memset(cross, 0, WORKER_NUM * TASK_NUM * sizeof(float));
		memset(fnet, 0, NN*NN * sizeof(uint8_t));
		memset(adj, 0, NN*NN * sizeof(short int));
		memset(deg, 0, NN * sizeof(short int));
		memset(par, 0, NN * sizeof(short int));
		memset(q, 0, NN * sizeof(short int));
		memset(inq, 0, NN * sizeof(short int));
		qs = 0;
		memset(d, 0, NN * sizeof(float));
		memset(pi, 0, NN * sizeof(float));
	}
};
