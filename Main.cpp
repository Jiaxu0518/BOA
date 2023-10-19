// Hybride.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "Opt.h"
#include "Greedy_RT.h"
#include "Greedy.h"
#include "Greedy_OT.h"
#include "Greedy_OTVar.h"
#include "SegOpt.h"
#include "SegFixOpt.h"
#include "SegVarOpt.h"
#include "SegBudOpt.h"
#include "SegBudFixOpt.h"
#include "SegBudVarOpt.h"
#include "OTSegBud.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#endif


// The one and only application object

CWinApp theApp;

using namespace std;

int main(int argc, char* argv[])
{
    int nRetCode = 0;

    HMODULE hModule = ::GetModuleHandle(nullptr);

    if (hModule != nullptr)
    {
        // initialize MFC and print and error on failure
        if (!AfxWinInit(hModule, nullptr, ::GetCommandLine(), 0))
        {
            // TODO: change error code to suit your needs
            wprintf(L"Fatal Error: MFC initialization failed\n");
            nRetCode = 1;
        }
		else
		{
			// TODO: code your application's behavior here.
			//为全局变量赋值			
			mode = atoi(argv[1]); //true:syn false real
			path.assign(argv[2]); //数据文件路径
			wfile.assign(path + argv[3]); //工人数据文件
			tfile.assign(path + argv[4]); //任务数据文件
			ofile.assign(path + argv[5]); //offline最优匹配结果文件
			datfile.assign(path + argv[6]); //结果统计文件
			sheet = atoi(argv[7]);  //数据保存到sheet 编号 = 1; 
			al = atoi(argv[8]);     //运行的算法编号 =	1
			B = atof(argv[9]); //总预算 syn:3000 real:300
			velocity = atof(argv[10]); //工人速度(km/h) syn:1 real:40
			cmax = atof(argv[11]);	//工人到任务最大距离 syn:1000 real:41.7027distance between (40.5998,-74.0701) and (40.8998,-73.7701)
			SEGS = atoi(argv[12]);   //时间段数，syn:10 real:72段,每段10分钟
			TIME = atoi(argv[13]);  //整个timeline长度，syn:100 real:720,表示720分钟
			
			CoInitialize(NULL); //EXCEL初始化
			switch (al)
			{
				case 1:/******************Opt*************************/
				{
					//检查opt文件是否存在，如果不存在创建一个空文件，如果存在的话清空。
					Util::CreateExcelFile(ofile.c_str());
					Info info_opt;
					Opt *opt = new Opt(B);
					opt->run(info_opt);
					delete opt;
					Util::save_statistic(datfile.c_str(), sheet, info_opt);
					break;
				}
				case 2:/*****************Greedy-RT Single**************/
				{					
					Info info_ng;
					double threshold = exp(Util::unirand(0, (int)log(cmax + 1)));
					cout << "cmax:" << cmax << " threshold:" << threshold << endl;
					Greedy_RT rt(B, threshold);
					rt.run_single(info_ng);
					Util::save_statistic(datfile.c_str(), sheet, info_ng);
					break;
				}
				case 3:/*****************Greedy-RT Expectation**************/
				{	
					Info info_ng;
					double count = 0;						
					double cost = 0;
					double totaltime = 0;
					double totoalmem = 0;
					int n = (int)log(cmax + 1);
					for (int i = 0; i <=n; i++)
					{	
						double threshold = exp(i);
						cout << "cmax:" << cmax << " threshold:" << threshold <<endl;
						Greedy_RT rt(B, threshold);
						vector<pair<pair<int, int>, double>> mat = rt.run_ave(totaltime,totoalmem);
						count += mat.size();
						for (vector<pair<pair<int, int>, double>>::iterator it = mat.begin();
							it != mat.end(); it++)
							cost += (*it).second;
					}
					Util::setInfo(info_ng, "Greedy-RT_AVE", cost/(n+1), totaltime/(n+1), totoalmem/(n+1), count/(n+1));
					Util::save_statistic(datfile.c_str(), sheet, info_ng);
					break;
				}
				case 4:/*************Greedy Algorithm***********************/
				{					
					Info info_greedy;
					Greedy g(B);
					g.run(info_greedy);
					Util::save_statistic(datfile.c_str(), sheet, info_greedy);
					break;
				}
				case 5:/*************Greedy_OT Algorithm***********************/
				{				
					Info info_greedy_ot;
					Util::loadOpt(ofile.c_str());
					double threshold = Util::getOptMaxDis(opt);
					cout << "OPT Threshold:" << threshold << endl;
					Greedy_OT got(B,threshold);
					got.run(info_greedy_ot);
					Util::save_statistic(datfile.c_str(), sheet, info_greedy_ot);
					break;
				}
				case 6:/*************Greedy_OTVar Algorithm*******各时间段采用的阈值不同****************/
				{					
					//动态创建数组存放各时间段的阈值和匹配数量
					double *threshold = new double[SEGS];
					memset(threshold, 0, SEGS * sizeof(double));
					int *count = new int[SEGS];
					memset(count, 0, SEGS * sizeof(int));
					Util::getOptSegMaxDis(ofile.c_str(), threshold, count, SEGS, TIME);

					Info info_greedy_ot;
					Greedy_OTVar gotvar(B, threshold,SEGS,TIME);
					gotvar.test(info_greedy_ot);
					Util::save_statistic(datfile.c_str(), sheet, info_greedy_ot);

					delete[] threshold;
					delete[] count;
					break;
				}
				case 7:/******************Segs_Opt********每个时间段做一次最优匹配，性能正常应该比较差*****************/
				{
					Info info_opt;
					SegOpt *segopt = new SegOpt(B,SEGS,TIME);
					segopt->test(info_opt);
					Util::save_statistic(datfile.c_str(), sheet, info_opt);
					delete segopt;
					break;
				}
				case 8:/******************Segs_Fix_Opt******在每个时间段根据最大阈值筛选掉travel cost大的工人后，再做最优匹配*******************/
				{				
					Info info_opt;
					//获取最优阈值
					Util::loadOpt(ofile.c_str());
					double threshold = Util::getOptMaxDis(opt);

					SegFixOpt *segfixopt = new SegFixOpt(B, SEGS, TIME,threshold);
					segfixopt->test(info_opt);
					Util::save_statistic(datfile.c_str(), sheet, info_opt);
					delete segfixopt;
					break;
				}
				case 9:/******************Segs_Var_Opt*********在每个时间段根据该时间段最大阈值筛选工人后，再做最优匹配****************/
				{					
					Info info_opt;
					SegVarOpt *segvaropt = new SegVarOpt(B, SEGS, TIME);
					segvaropt->test(info_opt);
					Util::save_statistic(datfile.c_str(), sheet, info_opt);
					delete segvaropt;
					break;
				}
				case 10:/******************Segs_Bud_Opt******每个时间段做最优匹配（各时间段投入预算不同）*******************/
				{					
					Info info_opt;
					SegBudOpt *segbudopt = new SegBudOpt(B, SEGS, TIME);
					segbudopt->test(info_opt);
					Util::save_statistic(datfile.c_str(), sheet, info_opt);
					delete segbudopt;
					break;
				}
				case 11:/******************Segs_Bud_Fix_Opt*****分段预算，各时间段先用最大阈值筛选完工人后再做最优匹配**********/
				{					
					Info info_opt;
					//获取最优阈值
					Util::loadOpt(ofile.c_str());
					double threshold = Util::getOptMaxDis(opt);

					SegBudFixOpt *sbfo = new SegBudFixOpt(B, SEGS, TIME, threshold);
					sbfo->test(info_opt);
					Util::save_statistic(datfile.c_str(), sheet, info_opt);
					delete sbfo;
					break;
				}
				case 12:/******************Segs_Bud_Var_Opt*******分段预算，各时间段先用分段最大阈值筛选完工人后再做最优匹配******************/
				{
					
					Info info_opt;
					SegBudVarOpt *sbvo = new SegBudVarOpt(B, SEGS, TIME);
					sbvo->test(info_opt);
					Util::save_statistic(datfile.c_str(), sheet, info_opt);
					delete sbvo;
					break;
				}
				case 13:/******************OT_Segs_Bud*****分段预算，采用单一最优阈值贪心算法完成任务配**********/
				{
					Info info_opt;
					//获取最优阈值
					Util::loadOpt(ofile.c_str());
					double threshold = Util::getOptMaxDis(opt);

					OTSegBud osb(B, SEGS, TIME, threshold);
					osb.test(info_opt);
					Util::save_statistic(datfile.c_str(), sheet, info_opt);
					break;
				}
				case 14: 
				{	/*************Generate Adversary workers***********************/
					//Util::CreateExcelFile(disfile.c_str());
					//Util::generateAdversaryWorkers(wfile.c_str(),tfile.c_str(),disfile.c_str());
					//Util::getMinDisForAllWorkers(wfile.c_str(), tfile.c_str(), disfile.c_str());
					/*******************计算两点之间的gps距离***********************/
					//worker[0].pos_x = -74.0701; task[0].pos_x = -73.7734;
					//worker[0].pos_y = 40.5998; task[0].pos_y = 40.8995;
					//worker[0].appear = 0; worker[0].id = 1;
					//task[0].id = 1; task[0].release = 0; task[0].deadline = 10000;
					//cout << "dis " << Util::Euclidean_Real(1, 1) << endl;
					//cout << "distance:"<<Util::Manhattan_Real(-73.9201, 40.8998 ,-73.9201, 40.5998) << endl;
					/******************测试各时间段最大阈值***********************/
					//double p[72] = {};//数组初始化0
					//int count[72] = {};
					//double bud[72] = {};
					//Util::getOptSegMaxDis(ofile.c_str(),p,count,72,720);
					//for (int i = 0; i < 72; i++)
					//	printf("第 %d 个时间段最大阈值 %f,数量%d\n", i,p[i],count[i]);
					//Util::getOptSegBudget(ofile.c_str(), bud, 72, 720);
					//for (int i = 0; i < 72; i++)
					//	printf("第 %d 个时间段预算 %f\n", i, bud[i]);
					/***************测试Rtree***********************************/
					/*Util::loadWorkers(wfile.c_str());
					Util::loadTasks(tfile.c_str());
					cout << "构建任务Rtree中..."<<endl;
					Util::buildRtree();
					cout << "查找最近邻..." << endl;
					for (int wid = 1; wid <= WORKER_NUM; wid++)
					{
						pair<pair<int, int>, double> p = Util::nearestTask_rt(wid);
						cout << "worker" << p.first.first << "'s nearest taskid:" << p.first.second << " distance:" << p.second << endl;
					}*/
					//for (int wid = 1; wid <= WORKER_NUM; wid++)
					//{
					//	double mind = FINF;
					//	int mtid = -1;
					//	for (int i = 1; i <= TASK_NUM; i++)
					//	{
					//		if (fabs(Util::Euclidean_Real(wid, i) + 1) <= eps) continue;
					//		if (Util::greater(mind, Util::Euclidean_Real(wid, i), 0))
					//		{
					//			mind = Util::Euclidean_Real(wid, i);
					//			mtid = i;
					//		}
					//	}
					//	cout << "worker" << wid << "'s nearest taskid:" << mtid << " distance:" << mind << endl;
					//}
					
					//cout << "输出Rtree..." << endl;
					////Util::OutputRTree();
					break;
				}
			}
		}
    }
    else
    {
        // TODO: change error code to suit your needs
        wprintf(L"Fatal Error: GetModuleHandle failed\n");
        nRetCode = 1;
    }

    return nRetCode;
}