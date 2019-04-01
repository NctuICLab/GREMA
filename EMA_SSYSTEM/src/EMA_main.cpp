//
//  EMA_main.cpp
//  EMA_ssystem
//
//  Created by Mingju on 2015/6/3.
//  Copyright (c) 2015年 Mingju. All rights reserved.
//
#include "stdafx.h"
#include "IGA.h"
#include "OrthogonalArray.h"
#include "S-system-div.h"
#include "getknowledge.h"
#include "Objfunction.h"
using namespace std;
int iteration;
int test_mode,cc_mode;
int ExecutionRuns;
int FUNC;
int GeneIndex;
int POP_SIZE;
int GENERATION;
int NO_timepoints;
int init_file_idx;
int g_knowledge_idx;
int h_knowledge_idx;
//test
void exit_with_help();
void parse_command_line(int argc, char **argv);

int main(int argc, char **argv)  
{
    Knowledge *knowledge;
    S_System_Divide *s_system;
    COrthogonalArray *OA;
    IGA *iga;
    Objfunction *obj;
    parse_command_line(argc, argv);
    unsigned int seed;
    seed = (unsigned int)time(NULL);
    srand(seed);
    //char DataFileName[30];
    int i,j;
    int time_step;
    unsigned int execution_times;
    unsigned int NumGeneParam;
	char str1[100],str2[100];
    FILE *ptr1,*ptr2,*ptr3;
    sprintf(str1,"IGA_result%d_step%d.txt",GeneIndex,iteration);
    ptr2=fopen(str1,"w");
    sprintf(str2,"calprofile%d_step%d.txt",GeneIndex,iteration);
    ptr3=fopen(str2,"w");
    
    s_system = new S_System_Divide();
    //strcpy(DataFileName,"./data/GNdata.txt");
    s_system->ReadInitFromFile(argv[init_file_idx]);
    s_system->CalAllExpValues(); //produce simulated data and perturbation set
    s_system->SetGeneIndex(GeneIndex);
    ptr1=fopen("printdata.txt","w");
    s_system->PrintAllExpXs(ptr1);//print expression data
    if(s_system->IsSimulated()){time_step = s_system->NumTrials()/NO_timepoints;}
    else{time_step = 1;}
    //============
    knowledge = new Knowledge();
    knowledge->NumVars = s_system->NumVariables();
    knowledge->NumTFs = s_system->NumTFsVariables();
    knowledge->initial_knowledge();
    knowledge->Readknowledge(GeneIndex,argv[g_knowledge_idx],argv[h_knowledge_idx]);
    //============
    obj = new Objfunction(s_system,knowledge);
    obj->NumVars = s_system->NumVariables();
    obj->NumTFs = s_system->NumTFsVariables();
    obj->NumTrials = s_system->NumTrials();
    obj->NumRuns = s_system->NumRuns();
    //============
    OA = new COrthogonalArray(obj);
    iga = new IGA(OA,s_system,knowledge,obj);
    iga->NumVars = s_system->NumVariables();
    iga->NumTFs = s_system->NumTFsVariables();
    iga->NumRuns = s_system->NumRuns();
    iga->NumTrials = s_system->NumTrials();
    iga->Num_nochange = s_system->NumNochange();
    iga->NO_g = knowledge->get_gNumconnections();
    iga->NO_h = knowledge->get_hNumconnections();
    iga->NO_POP = POP_SIZE;
    iga->NO_GEN = GENERATION;
    iga->set_parameter(GeneIndex,iteration,ExecutionRuns);
    iga->initialIGA();
    iga->RangeVariable();
    //============
    obj->PARAM_NUM = iga->get_PARAM_NUM();
    obj->NumGeneParam =iga->get_NumGeneParam();
    obj->NO_POP = POP_SIZE;
    obj->init_obj(time_step,FUNC,ExecutionRuns);

    NumGeneParam = iga->get_NumGeneParam();
    
    
//----------------------------------------------------------------------------
    fprintf(ptr3,"Seed:%d\n",seed);
    for(execution_times=0;execution_times<ExecutionRuns;execution_times++)
    {
        iga->run(execution_times);
        for(i=0;i<s_system->NumRuns();i++)
        {
            fprintf(ptr3,"N:%d Run:%d ",execution_times,i);
            for(j = 0; j < s_system->NumTrials(); j += time_step)
            {
                fprintf(ptr3,"%g ",obj->cal_profile[i][j]);
            }//end of NumTrials
            fprintf(ptr3,"%e %1.3f\n",obj->EVAL[iga->get_BestOne()],obj->best_cc);
        }
    }
    for(i=0;i<ExecutionRuns;i++)
    {
        for(j=0;j<=NumGeneParam;j++)
        {
            if(j==NumGeneParam)
                fprintf(ptr2,"%e ",obj->ModelParameters[i][j]);
            else
                fprintf(ptr2,"%f ",obj->ModelParameters[i][j]);
        }
        fprintf(ptr2,"%1.3f\n",obj->best_cc);
    }
    fclose(ptr1);
    fclose(ptr2);
    fclose(ptr3);
    return 0;
}

void exit_with_help()
{
	printf(
           "Usage: EMA-ssystem [options] init_file gknowledge_file hknowledge_file\n"
           "options:\n"
           "-i : set no of gene\n"
           "-I : set no of Iteration\n"
           "-n : set number of EMA's run (default 30)\n"
           "-G : set number of generation in EMA (default 10000)\n"
           "-F : set objective function in EMA (default 2)\n"
           "    0 -- (Xcal-Xexp)/Xexp\n"
           "    1 -- (Xcal-Xexp)\n"
           "    2 -- delta_diff + f0\n"
           "    3 -- delta_diff + f1\n"
           "    4 -- delta_diff\n"
           "-P : set number of popluation size (default 100)\n"
           "-t : set time points (necessary in simulated exp., default 10)\n"
           "-m change_mode : set the test mode or run mode (default 1)\n"
           "	0 -- run mode\n"
           "	1 -- test mode\n"
           "-c cc_mode : set the fitness function with cc or not (default 0)\n"
           "	0 -- without cc mode\n"
           "	1 -- with cc mode\n"           );
	exit(1);
}
void parse_command_line(int argc, char **argv)
{
	int i;
    	if(argc <= 1)
        	exit_with_help();
	//===default values=====//
	GeneIndex = 0;
	iteration = 1;
	ExecutionRuns = 30;
	GENERATION = 10000;
	NO_timepoints = 10;
	POP_SIZE = 100;
	FUNC = 2;
	test_mode = 1;
	cc_mode = 0;
	init_file_idx = 2;
	g_knowledge_idx = 3;
	h_knowledge_idx = 4;
        //======================//
    	i = 1;
	while(1)
	{
		if(argv[i][0] != '-')
			break;
		switch(argv[i][1])
		{
			case 'i':
				GeneIndex = atoi(argv[i+1]);
				i++;
				break;
            		case 'I':
                		iteration = atoi(argv[i+1]);
				i++;
				break;
			case 'n':
				ExecutionRuns = atoi(argv[i+1]);
				i++;
				break;
			case 'G':
				GENERATION = atoi(argv[i+1]);
				i++;
				break;
			case 'F':
				FUNC = atoi(argv[i+1]);
				i++;
				break;
            		case 'P':
				POP_SIZE = atoi(argv[i+1]);
				i++;
				break;
			case 't':
				NO_timepoints = atoi(argv[i+1]);
				i++;
				break;
            		case 'm':
                		test_mode = atoi(argv[i+1]);
				i++;
				break;
            		case 'c':
				cc_mode = atoi(argv[i+1]);
				i++;
				break;
			default:
				fprintf(stderr,"Unknown option!!\n");
				exit_with_help();
				break;
		}
		i++;
	}
	if(i < argc)
		init_file_idx = i++;
	else{
		fprintf(stderr, "No initial file!!\n");
		exit_with_help();
	}
	if(i < argc)
		g_knowledge_idx = i++;
	else{
		fprintf(stderr, "No g knowledge file!!\n");
		exit_with_help();
	}
	if(i < argc)
		h_knowledge_idx = i;
	else{
		fprintf(stderr, "No h knowledge file!!\n");
		exit_with_help();
	}
}
