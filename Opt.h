#pragma once
#include "Util.h"
#include "psapi.h"
#define CLR(a, x) memset( a, x, sizeof( a ) )
#define BUBL { \
    t = q[i]; q[i] = q[j]; q[j] = t; \
    t = inq[q[i]]; inq[q[i]] = inq[q[j]]; inq[q[j]] = t; }

#define Pot(u,v) (d[u] + pi[u] - pi[v])
#pragma comment(lib, "Psapi.lib")
using namespace std;
class Opt {
private:
	double budget = 0;
	uint8_t cap[NN][NN];
	// cost per unit of flow matrix (fill this up)
	float s[WORKER_NUM];
	float e[TASK_NUM];
	float cross[WORKER_NUM][TASK_NUM];
	// flow network and adjacency list
	uint8_t fnet[NN][NN];
	short int adj[NN][NN], deg[NN];
	// Dijkstra's predecessor, depth and priority queue
	short int par[NN], q[NN], inq[NN], qs;
	float d[NN];
	float pi[NN];
public:
	Opt(double budget)
	{
		this->budget = budget;
	}
	~Opt()
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
	void run(Info &item)
	{
		//加载工人和任务数据
		Util::loadWorkers(wfile.c_str());
		Util::loadTasks(tfile.c_str());
		//图中节点数,包括源和终节点
		int numV = WORKER_NUM + TASK_NUM + 2;
		memset(cap, 0, sizeof(cap));//clear 0
		int s = 0, t = numV - 1;
		//构建图
		buildGraph(WORKER_NUM, TASK_NUM);
		float fcost;
		PROCESS_MEMORY_COUNTERS cos;
		clock_t start, finish;
		cout << "matching...." << endl;
		start = clock();
		int flow = mcmf4(numV, s, t, fcost, budget);
		finish = clock();
		::GetProcessMemoryInfo(GetCurrentProcess(), &cos, sizeof(cos));
		//提取并保存匹配结果
		vector<pair<pair<int, int>, double>> mat = optMatching();
		Util::save(mat, ofile.c_str());
		Util::setInfo(item, "OPT", (double)fcost, (double)(finish - start) / CLOCKS_PER_SEC, cos.WorkingSetSize / 1024, mat.size());
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
		fcost = 0;
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
	//搭建图
	void buildGraph(int workerNum, int taskNum)
	{
		cout << "building grpah..." << endl;
		for (int i = 1; i <= workerNum; i++) {
			s[i - 1] = 0; //费用
			cap[0][i] = 1;  //容量
		}
		for (int i = workerNum + 1; i <= workerNum + taskNum; i++) {
			e[i - workerNum - 1] = 0; 
			cap[i][workerNum + taskNum + 1] = 1; 
		}
		for (int taskid = workerNum + 1; taskid <= workerNum + taskNum; taskid++) {
			for (int workerid = 1; workerid <= workerNum; workerid++)
			{
				int wid = workerid;
				int tid = taskid - workerNum;
				double dis;
				if (mode)
					dis = Util::Manhattan(wid, tid);  //人工数据集
				else
					dis = Util::Euclidean_Real(wid, tid); //真实数据集
				if (fabs(dis + 1)>eps) //工人在dealine前赶到任务地点
				{
					cross[workerid - 1][taskid - workerNum - 1] = dis; // cross数组存放的是工人和任务之间的cost;
					cap[workerid][taskid] = 1; 
				}
				else
				{
					cross[workerid - 1][taskid - workerNum - 1] = FINF;
				}
			}
		}
	}
	//提取匹配结果
	vector<pair<pair<int, int>, double>> optMatching()
	{
		vector<pair<pair<int, int>, double>> v;
		for (int windex = 1; windex <= WORKER_NUM; windex++)
			for (int tindex = WORKER_NUM + 1; tindex <= WORKER_NUM + TASK_NUM; tindex++)
			{
				if (fnet[windex][tindex] == 1)
					v.push_back(make_pair(make_pair(windex, tindex), cost(windex, tindex))); //保存工人和任务的图ID
			}
		return v;
	}
};
