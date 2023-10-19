#pragma once
//#include <vector>
#include <iostream>
#include <IllusionExcelFile.h>
#include <time.h>
#include <cstdlib>
#include <random>
//#include <algorithm>
//#include <cmath>
#include "rtree.h"
#define FINF 10000
#define eps 1e-5
#define EARTH_RADIUS 6378.137
#define M_PI 3.14159265358979323846
using namespace std;
using namespace rbox;

struct Worker  //id，x-coordinate，y-coordinate，appear time，deadline
{
	int id;   
	double pos_x;
	double pos_y;
	time_t appear;    
	time_t deadline;  
};
struct Task //id, x-coordinate, y-coordinate, release time, deadline
{
	int id;
	double pos_x;
	double pos_y;
	time_t release;  
	time_t deadline; 
	int matched = 0;
};
struct Info
{
	char *a_name;
	double cost; //used budget
	double time; //run time
	double memory; //used memory
	int matched;  //matching size
	int v_worker; //number of workers
	int v_task; //number of tasks
	double v_budget; //total budget
	int nsegs;
};
//rtee矩形区域
struct Rect
{
	Rect() {}

	Rect(double a_minX, double a_minY, double a_maxX, double a_maxY)
	{
		min[0] = a_minX;
		min[1] = a_minY;

		max[0] = a_maxX;
		max[1] = a_maxY;
	}


	double min[2];
	double max[2];
};
// the number of workers and tasks
#define NN 5682+6000+2  //设置最优算法节点数目
static const int WORKER_NUM = 5682;
static const int TASK_NUM = 6000;
static Worker worker[WORKER_NUM];  //can not swap or sort the position in array
static Task task[TASK_NUM];  //can not change the position in array
string path;//数据路径 = "E:\\temp\\experiment\\data\\synthetic_final"
string wfile; //工人数据文件 = path + "\\worker6000.xlsx"
string tfile; //任务数据文件 = path + "\\task6000.xlsx"
string ofile; //最优匹配方案存放位置 = path + "\\opt_t6000.xlsx"
string datfile; //结果汇总文件 = path + "\\dattest.xlsx"
string disfile; //测试工人和任务最近pair文件 = path + "\\testdis.xlsx"

 bool mode; //true:syn false real
 int sheet=-1;  //数据保存到sheet 编号 = 1; 
 int al=-1;     //运行的算法编号 =	1
 double B=-1; //total budget syn:3000 real:300
 double velocity = -1; //(km/h) syn:1 real:40
 double cmax = -1;	//syn distance syn:1000 real:41.7027distance between (40.5998,-74.0701) and (40.8998,-73.7701)
 int SEGS = -1;   //分段数，syn:10 real:72段,每段10分钟
 int TIME = -1;  //整个timeline长度，syn:100 real:720,表示720分钟
static RTree<int,double,2,double> tree; //int是值类型，double是坐标类型，2表示2维空间
static vector<pair<pair<int, int>, double>> opt;

class Util
{
public:

public:
	Util();
	~Util();

	static void loadWorkers(const char *wfile);
	static void loadTasks(const char *tfile);
	static void loadOpt(const char *optfile);
	static void display(vector<pair<pair<int, int>, double>> &v);
	static void save(vector<pair<pair<int, int>, double>> &mat, const char *file);
	static double Manhattan(int workerId, int taskId);
	static bool greater(double x, double y, int mode);//mode=1:>=,0:>
	static int unirand(int lb, int rb);
	static double getOptMaxDis(vector<pair<pair<int, int>, double>> &opt);
	static void save_statistic(const char *file, int sheet, Info item);
	static void setInfo(Info &item, char *algor_name, double cost, double time, double memory, int matched);
	static void setInfo(Info &item, char *algor_name, double cost, double time, double memory, int matched,int nsegs);
	static bool GreaterSort(pair<pair<int, int>, double> a, pair<pair<int, int>, double> b);
	static void generateAdversaryWorkers(const char *wfile, const char *tfile, const char *disfile);
	static void getMinDisForAllWorkers(const char *wfile, const char *tfile, const char *disfile);
	static void CreateExcelFile(const char *optfile);
	static double rad(double d);
	static double Euclidean_Real(int wid, int tid);
	static double Euclidean_Real(double min_x, double min_y, double max_x, double max_y);
	//获取各个时间段的最大阈值和匹配数量
	static void getOptSegMaxDis(const char *optfile,double thres[],int count[], int n,int time);
	//获取各个时间段分配的预算
	static void getOptSegBudget(const char *optfile, double bud[], int n, int time,double totalbudget);
	//根据task[TASK_NUM]构建任务Rtree
	static void buildRtree();
	//输出Rtree
	static void OutputRTree();
	//在RTree中寻找最近任务
	static pair<pair<int, int>, double> nearestTask_rt(int wid);
	static pair<pair<int, int>, double> nearestTask_rt(int wid, double threshold);
};

Util::Util()
{
}

Util::~Util()
{
}

double Util::rad(double d)
{
	return d * M_PI / 180.0;
}

void Util::loadWorkers(const char *wfile)
{
	//load workers from .xlsx file
	printf("Load workers...\n");
	IllusionExcelFile ief;
	if (ief.InitExcel())
	{ //init excel server
		if (ief.OpenExcelFile(wfile))
		{//open excel file
			ief.LoadSheet(1, true);//load sheet1
			for (int i = 1; i <= WORKER_NUM; i++) //construct RTree
			{
				worker[i - 1].id = ief.GetCellInt(i, 1);
				worker[i - 1].pos_x = ief.GetCellDouble(i, 2);
				worker[i - 1].pos_y = ief.GetCellDouble(i, 3);
				worker[i - 1].appear = ief.GetCellInt(i, 4);
				worker[i - 1].deadline = ief.GetCellInt(i, 5);
			}
		}
		else
			printf("open this file failed!\n");
		ief.CloseExcelFile(false);// close excel file
	}
	else
		printf("excell has not installed!\n");
	ief.ReleaseExcel();
}
void Util::loadTasks(const char *tfile)
{
	//load tasks from .xlsx file
	printf("Load tasks...\n");
	IllusionExcelFile ief;
	if (ief.InitExcel())
	{ //init excel server
		if (ief.OpenExcelFile(tfile))
		{//open excel file
			ief.LoadSheet(1, true);//load sheet1
			for (int i = 1; i <= TASK_NUM; i++) //construct RTree
			{
				task[i - 1].id = ief.GetCellInt(i, 1);
				task[i - 1].pos_x = ief.GetCellDouble(i, 2);
				task[i - 1].pos_y = ief.GetCellDouble(i, 3);
				task[i - 1].release = ief.GetCellInt(i, 4);
				task[i - 1].deadline = ief.GetCellInt(i, 5);
				task[i - 1].matched = 0;
			}
		}
		else
			printf("open this file failed!\n");
		ief.CloseExcelFile(false);// close excel file
	}
	else
		printf("excell has not installed!\n");
	ief.ReleaseExcel();
}
//load opt
void Util::loadOpt(const char *optfile)
{
	//load workers from .xlsx file
	printf("Load opt...\n");
	int optcount = 0;
	IllusionExcelFile ief;
	if (ief.InitExcel())
	{ //init excel server
		if (ief.OpenExcelFile(optfile))
		{//open excel file
			ief.LoadSheet(1, true);//load sheet1
			optcount = ief.GetRowCount();
			for (int i = 1; i <= optcount; i++) 
			{
				pair<pair<int, int>, double> p;
				p.first.first = ief.GetCellInt(i, 1);
				p.first.second = ief.GetCellInt(i, 2);
				p.second = ief.GetCellDouble(i, 3); //travel cost 
				opt.push_back(p);
			}
		}
		else
			printf("open this file failed!\n");
		ief.CloseExcelFile(false);// close excel file
	}
	else
		printf("excell has not installed!\n");
	ief.ReleaseExcel();
}
void Util::CreateExcelFile(const char *optfile)
{
	//load workers from .xlsx file
	printf("Create new File...\n");
	IllusionExcelFile ief;
	if (ief.InitExcel())
	{ //init excel server
		if (ief.OpenEmptyExcelFile(optfile)) //打开一个空文件
		{
			CString str(optfile);
			ief.SaveasXSLFile(str); //保存文件
		}
		else
			printf("open this file failed!\n");
		ief.CloseExcelFile(false);// close excel file
	}
	else
		printf("excell has not installed!\n");
	ief.ReleaseExcel();
}
//display optimal matching pair
void Util::display(vector<pair<pair<int, int>, double>> &v)
{
	vector<pair<pair<int, int>, double>>::iterator it;
	for (it = v.begin(); it != v.end(); it++)
	{
		int wid = (*it).first.first;
		int tid = (*it).first.second;
		cout << "w" << wid << "->" << "t" << tid << " cost:" << (*it).second << endl;
	}
}
//save matching result as opt_XXXXX.xlsx.
void Util::save(vector<pair<pair<int, int>, double>> &mat, const char *file)
{
	IllusionExcelFile ief;
	if (ief.InitExcel())
	{ //init excel server
		if (ief.OpenExcelFile(file))
		{
			ief.LoadSheet(1, true);//load sheet1,
			long row = 1;
			for (vector<pair<pair<int, int>, double>>::iterator it = mat.begin();
				it != mat.end(); it++, row++)
			{
				//transform id in graph to logical id
				ief.SetCellInt(row, 1, (*it).first.first);
				ief.SetCellInt(row, 2, (*it).first.second - WORKER_NUM);
				ief.SetCellDouble(row, 3, (*it).second);
				ief.SetCellInt(row, 4, worker[(*it).first.first-1].appear); //加入配对发生时间信息
			}
			CString str(file);
			ief.SaveasXSLFile(str);
		}
		else
			printf("open this file failed!\n");
		ief.CloseExcelFile(false);// close excel file
	}
	else
		printf("excell has not installed!\n");
	ief.ReleaseExcel();
}
void Util::save_statistic(const char *file, int sheet, Info item)
{
	IllusionExcelFile ief;
	if (ief.InitExcel())
	{ //init excel server
		if (ief.OpenExcelFile(file))
		{
			ief.LoadSheet(sheet, true);//load sheet,
			long row = ief.GetRowCount() + 1; //if row number is 0, return 1,this is the bug of getRowCount function.
			CString algo(item.a_name);
			ief.SetCellString(row, 1, algo);
			ief.SetCellDouble(row, 2, item.cost);
			ief.SetCellDouble(row, 3, item.time);
			ief.SetCellDouble(row, 4, item.memory);
			ief.SetCellInt(row, 5, item.matched);
			ief.SetCellInt(row, 6, item.v_worker);
			ief.SetCellInt(row, 7, item.v_task);
			ief.SetCellDouble(row, 8, item.v_budget);
			ief.SetCellInt(row, 9, item.nsegs);
			CString str(file);
			ief.SaveasXSLFile(str);
		}
		else
			printf("open this file failed!\n");
		ief.CloseExcelFile(false);// close excel file
	}
	else
		printf("excell has not installed!\n");
	ief.ReleaseExcel();
}
void Util::setInfo(Info &item, char *algor_name, double cost, double time, double memory, int matched)
{
	item.a_name = algor_name;
	item.cost = cost;
	item.time = time;
	item.memory = memory;
	item.matched = matched;
	item.v_worker = WORKER_NUM;
	item.v_task = TASK_NUM;
	item.v_budget = B;
	item.nsegs = 1;
}
void Util::setInfo(Info &item, char *algor_name, double cost, double time, double memory, int matched,int nsegs)
{
	item.a_name = algor_name;
	item.cost = cost;
	item.time = time;
	item.memory = memory;
	item.matched = matched;
	item.v_worker = WORKER_NUM;
	item.v_task = TASK_NUM;
	item.v_budget = B;
	item.nsegs = nsegs;
}
//Manhattan distance,workerId and taskId are logical id.
double Util::Manhattan(int workerId, int taskId)
{	//transform logical id to array id
	double dis = fabs(worker[workerId-1].pos_x - task[taskId-1].pos_x) + fabs(worker[workerId-1].pos_y - task[taskId-1].pos_y);
	//if worker arrives the target task before deadline,return the true travel cost.Otherwise, return inf.
	if (greater(task[taskId - 1].deadline,worker[workerId - 1].appear + dis,1))
		return dis;
	else
		return -1; //worker can't arrive task before it's deadline.
}
//GPS Euclidean distance (Km)
double Util::Euclidean_Real(int workerId,int taskId)
{
	double radLat1 = rad(worker[workerId - 1].pos_y);
	double radLat2 = rad(task[taskId - 1].pos_y);
	double a = radLat1 - radLat2;
	double b = rad(worker[workerId - 1].pos_x) - rad(task[taskId - 1].pos_x);
	double dis = 2 * asin(sqrt(pow(sin(a / 2), 2) +
		cos(radLat1)*cos(radLat2)*pow(sin(b / 2), 2)));
	dis = dis * EARTH_RADIUS;
	dis = round(dis * 10000) / 10000;
	double time = dis / velocity * 60; //change hour to minitues
	if (greater(task[taskId - 1].deadline, worker[workerId - 1].appear + time, 1))
		return dis;
	else
		return -1; //worker can't arrive task before it's deadline.
}
//GPS Manhattan distance (Km)
double Util::Euclidean_Real(double min_x, double min_y, double max_x, double max_y)
{
	double radLat1 = rad(min_y);
	double radLat2 = rad(max_y);
	double a = radLat1 - radLat2;
	double b = rad(min_x) - rad(max_x);
	double dis = 2 * asin(sqrt(pow(sin(a / 2), 2) +
		cos(radLat1)*cos(radLat2)*pow(sin(b / 2), 2)));
	dis = dis * EARTH_RADIUS;
	dis = round(dis * 10000) / 10000;
	return dis;
}
bool Util::greater(double x, double y, int mode)
{
	if (fabs(x - y) <= eps)
	{
		if (mode == 1) // >=
			return true;
		else              //>
			return false;
	}
	else if (x > y)
		return true;
	else
		return false;
}
int Util::unirand(int lb,int rb)
{
	random_device rd;
	// Use Mersenne twister engine to generate pseudo-random numbers.
	mt19937 engine(rd());
	// (Note that the range is [inclusive, inclusive].)
	uniform_int_distribution<int> dist(lb, rb);
	return dist(engine);
}
double Util::getOptMaxDis(vector<pair<pair<int, int>, double>> &opt)
{
	double maxDis = -1;
	for (vector<pair<pair<int, int>, double>>::iterator it = opt.begin(); it != opt.end(); it++)
		if (Util::greater((*it).second, maxDis, 0))
			maxDis = (*it).second;
	return maxDis;
}
void Util::getOptSegMaxDis(const char *optfile, double thres[], int count[], int n, int time)
{
	//load workers from .xlsx file
	printf("static Seg Max Travel Cost...\n");
	int optcount = 0;
	IllusionExcelFile ief;
	if (ief.InitExcel())
	{ //init excel server
		if (ief.OpenExcelFile(optfile))
		{//open excel file
			ief.LoadSheet(1, true);//load sheet1
			optcount = ief.GetRowCount();
			for (int i = 1; i <= optcount; i++)
			{
				double cost = ief.GetCellDouble(i, 3); //travel cost 
				int t = ief.GetCellInt(i, 4);//匹配发生时间
				int mappedIndex = t*n/time;  //计算该时间段对应的元素下标
				if (Util::greater(cost, thres[mappedIndex], 0))
					thres[mappedIndex] = cost;
				count[mappedIndex]++;  //计数
			}
		}
		else
			printf("open this file failed!\n");
		ief.CloseExcelFile(false);// close excel file
	}
	else
		printf("excell has not installed!\n");
	ief.ReleaseExcel();
}

void Util::getOptSegBudget(const char *optfile, double bud[], int n, int time, double totalbudget)
{
	//load workers from .xlsx file
	printf("static Seg budget to be allocated...\n");
	int optcount = 0;
	double budget_sum = 0;
	IllusionExcelFile ief;
	if (ief.InitExcel())
	{ //init excel server
		if (ief.OpenExcelFile(optfile))
		{//open excel file
			ief.LoadSheet(1, true);//load sheet1
			optcount = ief.GetRowCount();
			for (int i = 1; i <= optcount; i++)
			{
				double cost = ief.GetCellDouble(i, 3); //travel cost 
				int t = ief.GetCellInt(i, 4);//匹配发生时间
				int mappedIndex = t*n / time;  //计算该时间段对应的元素下标
				bud[mappedIndex] += cost; //统计各时间段travel cost之和
				budget_sum += cost; //计算travel cost总和
			}
			for (int i = 0; i < n; i++)
				bud[i] = totalbudget * bud[i] / budget_sum;
		}
		else
			printf("open this file failed!\n");
		ief.CloseExcelFile(false);// close excel file
	}
	else
		printf("excell has not installed!\n");
	ief.ReleaseExcel();
}

bool Util::GreaterSort(pair<pair<int, int>, double> a, pair<pair<int, int>, double> b)
{ 
	return (a.second>b.second); 
}
void Util::generateAdversaryWorkers(const char *wfile, const char *tfile, const char *disfile)
{
	loadWorkers(wfile);
	loadTasks(tfile);
	vector<pair<pair<int, int>, double>> result;
	for (int wid = 1; wid <= WORKER_NUM; wid++)
	{
		double mindis = FINF;
		int mintid = -1;
		for (int tid = 1; tid <= TASK_NUM; tid++)
		{
			double dis = fabs(worker[wid - 1].pos_x - task[tid - 1].pos_x) + fabs(worker[wid - 1].pos_y - task[tid - 1].pos_y);
			if (Util::greater(mindis, dis, 0))
			{
				mindis = dis;
				mintid = tid;
			}
		}
		result.push_back(make_pair(make_pair(wid, mintid), mindis));
	}
	//将result中的pair按照cost降序排列
	std::sort(result.begin(), result.end(), GreaterSort);
	//load workers from .xlsx file
	printf("Loader disfile...\n");
	IllusionExcelFile ief;
	if (ief.InitExcel())
	{ //init excel server
		if (ief.OpenExcelFile(disfile))
		{//open excel file
			ief.LoadSheet(1, true);//load sheet1
			long row = 1;
			for (vector<pair<pair<int, int>, double>>::iterator it = result.begin(); it != result.end(); it++)
			{
				pair<pair<int, int>, double> p = (*it);
				int wid = p.first.first;
				ief.SetCellInt(row, 1, row);
				ief.SetCellDouble(row, 2, worker[wid-1].pos_x);//替换原来行工人出现的坐标
				ief.SetCellDouble(row, 3, worker[wid - 1].pos_y);
				ief.SetCellInt(row, 4, worker[row - 1].appear); //保留原来行的出现时间和deadline
				ief.SetCellInt(row, 5, worker[row - 1].deadline);
				row++;
			}
			CString str(disfile);
			ief.SaveasXSLFile(str);
		}
		else
			printf("open this file failed!\n");
		ief.CloseExcelFile(false);// close excel file
	}
	else
		printf("excell has not installed!\n");
	ief.ReleaseExcel();
}
void Util::getMinDisForAllWorkers(const char *wfile, const char *tfile, const char *disfile)
{
	loadWorkers(wfile);
	loadTasks(tfile);
	vector<pair<pair<int, int>, double>> result;
	for (int wid = 1; wid <= WORKER_NUM; wid++)
	{
		double mindis = FINF;
		int mintid = -1;
		for (int tid = 1; tid <= TASK_NUM; tid++)
		{
			double dis = fabs(worker[wid - 1].pos_x - task[tid - 1].pos_x) + fabs(worker[wid - 1].pos_y - task[tid - 1].pos_y);
			if (Util::greater(mindis, dis, 0))
			{
				mindis = dis;
				mintid = tid;
			}
		}
		result.push_back(make_pair(make_pair(wid, mintid), mindis));
	}
	//将result中的pair按照cost降序排列
	std::sort(result.begin(), result.end(), GreaterSort);
	//load workers from .xlsx file
	printf("Loader disfile...\n");
	IllusionExcelFile ief;
	if (ief.InitExcel())
	{ //init excel server
		if (ief.OpenExcelFile(disfile))
		{//open excel file
			ief.LoadSheet(1, true);//load sheet1
			long row = 1;
			for (vector<pair<pair<int, int>, double>>::iterator it = result.begin(); it != result.end(); it++)
			{
				pair<pair<int, int>, double> p = (*it);
				ief.SetCellInt(row, 1, p.first.first);
				ief.SetCellInt(row, 2, p.first.second);
				ief.SetCellDouble(row, 3, p.second);
				row++;
			}
			CString str(disfile);
			ief.SaveasXSLFile(str);
		}
		else
			printf("open this file failed!\n");
		ief.CloseExcelFile(false);// close excel file
	}
	else
		printf("excell has not installed!\n");
	ief.ReleaseExcel();
}
//根据任务构建Rtree，当有工人到来通过Rtree寻找最近任务
void Util::buildRtree() {
	tree.RemoveAll();//清空原树
	printf("Build RTree...\n");
	for (int i = 1; i < TASK_NUM; i++)
	{
		Rect new_rect(task[i - 1].pos_x, task[i - 1].pos_y, task[i - 1].pos_x, task[i - 1].pos_y);
		tree.Insert(new_rect.min, new_rect.max, task[i - 1].id);//i is task id	
	}
}
//输出任务树
void Util::OutputRTree() {
	int itIndex = 1;
	RTree<int,double,2,double>::Iterator it;//iterator
	for (tree.GetFirst(it); !tree.IsNull(it); tree.GetNext(it))
	{
		int value = tree.GetAt(it);
		double boundsMin[2] = {};
		double boundsMax[2] = {};
		it.GetBounds(boundsMin, boundsMax);
		printf("it[%d] %d = (%f,%f,%f,%f)\n", itIndex++, value, boundsMin[0], boundsMin[1], boundsMax[0], boundsMax[1]);
	}
}
//在RT中查找最近任务
pair<pair<int, int>, double> Util::nearestTask_rt(int wid)
{
	int tid = -1;
	double mindis = FINF;
	vector<int> search_out;
	double coord[2] = { worker[wid - 1].pos_x, worker[wid-1].pos_y };
	tree.KNN(coord, 1, search_out,wid,mindis,FINF);
	if (search_out.size() == 1)
	{
		tid = search_out[0];
	}		
	return make_pair(make_pair(wid, tid), mindis);
}
//在RT中查找最近任务
pair<pair<int, int>, double> Util::nearestTask_rt(int wid,double threshold)
{
	int tid = -1;
	double mindis = FINF;
	vector<int> search_out;
	double coord[2] = { worker[wid - 1].pos_x, worker[wid - 1].pos_y };
	tree.KNN(coord, 1, search_out, wid, mindis,threshold);
	//tree.KNN(coord, 1, search_out, wid, mindis);
	if (search_out.size() == 1)
	{
		tid = search_out[0];
	}
	return make_pair(make_pair(wid, tid), mindis);
}


