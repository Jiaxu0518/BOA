BOA is an open source project for the paper "Budget-Aware Online Task Assignment Based on Historical Temporal-Spatial Records in Spatial Crowdsourcing". 
It is implemented using the C++ programming language and the MFC programming framework. 
The running of the program requires the following major environment Settings:
1. Install the MFC extension in Visual Studio.
https://developer.aliyun.com/article/1144923
2. Disable Intellisense.
https://blog.csdn.net/qq_43533553/article/details/123373431

When writing the paper, we found that the algorithm names in the program code is relatively arbitrary.
So in the paper, we changed these algorithm names to unified names which are easy to understand. 
This leads to the inconsistency of algorithm names in the paper and in the program. 
The following are the mappings of these algorithms.
---------------------------------------------------
Code file in BOA project	Algorithm in manuscript
---------------------------------------------------
SegOpt.h                  SEGOPT
Greedy_OT.h               OTGRY
Greedy_OTVar.h            SEGOTGRY
SegVarOpt.h               SEGOTOPT
SegBudOpt.h               SEGBUDOPT
MBSGA.h                   SEGBUDGRY
Greedy.h                  GREEDY
TBSGA.h                   SEGOTBUDGRY
SegBudVarOpt.h            SEGOTBUDOPT
---------------------------------------------------
