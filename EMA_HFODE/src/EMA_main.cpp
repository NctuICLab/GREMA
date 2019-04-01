//
//  EMA_main.cpp
//  EMA_Hillsum_mask
//
//  Created by Mingju on 2017/3/10.
//  Copyright (c) 2015 Mingju. All rights reserved.
//
#include "Define.h"
#include "IGA.h"
#include "OrthogonalArray.h"
#include "Model.h"
#include "Getknowledge.h"
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
int knowledge_idx;
//test
void exit_with_help();
void parse_command_line(int argc, char **argv);
int main(int argc, char **argv)  
{
    parse_command_line(argc, argv);
    Getknowledge *knowledge;
    Model *hill;
    COrthogonalArray *OA;
    IGA *iga;
    Objfunction *obj;
    double	**Chromosomes;
    Chromosomes = NULL;
    double **ModelParameters;
    ModelParameters = NULL;
    double *PopFitness;
    PopFitness = NULL;
    PopFitness = new double[POP_SIZE];
    memset(PopFitness, 0, POP_SIZE*sizeof(double));
    unsigned int seed;
    seed = (unsigned int)time(NULL);
    srand(seed);
    //char DataFileName[30];
    int i,j;
    int time_step;
    int bestIndex;
    int NumReg = 0;
    int NumVars;
    int NumTFs;
    int NumEncodeParam;
    int TimePoints;
    int ExpRepeat;
    int *range = NULL;
    int *choose = NULL;
    int *LocMask = NULL;
    //unsigned int ChromSize;
    unsigned int execution_times;
    unsigned int NumOrigParam;
    //char str1[100];
    char str2[100];
    //FILE *ptr1,*ptr2,*ptr3;
    FILE *ptr_Inputprofile;
    FILE *ptr_EMAcalProfile;
    //sprintf(str1,"%d_step%d.txt",GeneIndex,iteration);
    //ptr2=fopen(str1,"w");
    sprintf(str2,"calprofile%d_step%d.txt",GeneIndex,iteration);
    ptr_EMAcalProfile=fopen(str2,"w");
    
    hill = new Model();
    //strcpy(DataFileName,"./data/GNdata.txt");
    hill->SetGeneIndex(GeneIndex);
    hill->ReadInitFromFile(argv[init_file_idx]);
    if(hill->IsSimulated()){
        hill->CalAllExpValues(); //produce simulated data and perturbation set
    }
    ptr_Inputprofile=fopen("printdata.txt","w");
    hill->PrintAllExpXs(ptr_Inputprofile);//print expression data
    //exit(0);
    if(hill->IsSimulated()){time_step = hill->NumTrials()/NO_timepoints;}
    else{time_step = 1;}
    NumVars = hill->NumVariables();
    NumTFs = hill->NumTFsVariables();
    TimePoints = hill->NumTrials();
    ExpRepeat = hill->NumRuns();
    range = new int [NumVars];
    memset(range,0,NumVars*sizeof(int));
    choose = new int [NumVars];
    memset(choose, 0, NumVars*sizeof(int));
    LocMask = new int [NumVars];
    memset(LocMask, 0, NumVars*sizeof(int));
    /*
    printf("initial\n");
    for(int i=0; i<NumVars; i++){
        printf("LockMask[%d]=%d\trange[%d]=%d\tchoose[%d]=%d\n",i,LocMask[i],i,range[i],i,choose[i]);
    }
    PAUSE;
    */
    //============
    knowledge = new Getknowledge();
    knowledge->NumVars = NumVars;
    knowledge->NumTFs = NumTFs;
    //knowledge->initial_knowledge();
    knowledge->Readknowledge(GeneIndex, iteration, LocMask, range, choose, argv[knowledge_idx]);
    /*
    printf("read knowledge of regulations file\n");
    for(int i=0; i<NumVars; i++){
        printf("LockMask[%d]=%d\trange[%d]=%d\tchoose[%d]=%d\n",i,LocMask[i],i,range[i],i,choose[i]);
    }
    PAUSE;
     */
    NumReg = knowledge->get_Numconnections();
    NumEncodeParam = 3*NumReg + 4;//Location, Nij, kij + bi, betai, degi
    NumOrigParam = 2*NumVars +3;//Nij,kij + bi, betai, degi
    //ChromSize = POP_SIZE * NumEncodeParam;
    Chromosomes = new double* [POP_SIZE];
    for(i=0; i<POP_SIZE; i++){
        Chromosomes[i] = new double [NumEncodeParam];
        memset(Chromosomes[i], 0, NumEncodeParam*sizeof(double));
    }
    ModelParameters = new double* [ExecutionRuns];
    for(i=0; i<ExecutionRuns; i++){
        ModelParameters[i] = new double [NumOrigParam+1];
        memset(ModelParameters[i], 0, (NumOrigParam+1)*sizeof(double));
    }
    /*
    printf("initial the chromosome value\n");
    for(int i=0; i<POP_SIZE; i++){
        for(int j=0; j<NumEncodeParam; j++){
            printf("%f ",Chromosomes[i][j]);
        }
        printf("\n");
    }
    PAUSE;
     */
    obj = new Objfunction();
    obj->NumVars = NumVars;
    obj->NumTFs = NumTFs;
    obj->NumTrials = TimePoints;
    obj->NumRuns = ExpRepeat;
    obj->PARAM_NUM = NumEncodeParam;
    obj->NumGeneParam = NumOrigParam;
    obj->NO_POP = POP_SIZE;
    obj->init_obj(time_step,FUNC,ExecutionRuns,cc_mode, NumReg,GeneIndex,hill,Chromosomes);
    //==========================================
    OA = new COrthogonalArray();
    OA->initialOA(obj, Chromosomes, PopFitness);
    //===========================================
    //iga = new IGA(OA,hill,knowledge,obj);
    iga = new IGA();
    iga->NumVars = NumVars;
    iga->NumTFs = NumTFs;
    iga->NumRuns = ExpRepeat;
    iga->NumTrials = TimePoints;
    iga->Num_nochange = hill->NumNochange();
    iga->NO_REG = NumReg;
    iga->NO_POP = POP_SIZE;
    iga->NO_GEN = GENERATION;
    iga->PARAM_NUM = NumEncodeParam;
    iga->NumGeneParam = NumOrigParam;
    iga->set_parameter(GeneIndex,iteration,test_mode,OA,hill,knowledge,obj,Chromosomes,PopFitness,LocMask,choose,range, ModelParameters);
    iga->initialIGA();
    iga->RangeVariable();
    
    
    
//----------------------------------------------------------------------------
    fprintf(ptr_EMAcalProfile,"Seed:%d\n",seed);
    for(execution_times=0;execution_times<ExecutionRuns;execution_times++)
    {
        if(knowledge->get_noregulation())
        {
            for(i=0; i<POP_SIZE; i++){
                for(j=0; j<NumEncodeParam; j++){
                    Chromosomes[i][j] = 0;
                    //printf("%f ",Chromosomes[i][j]);
                }
                //printf("\n");
            }
            
        }else{
            iga->initialpop();
            iga->run(execution_times);
        }
        for(i=0; i<ExpRepeat; i++)
        {
            fprintf(ptr_EMAcalProfile,"N:%d Run:%d ",execution_times,i);
            for(j=0; j<TimePoints; j+=time_step)
            {
                fprintf(ptr_EMAcalProfile,"%g ",obj->cal_profile[i][j]);
            }//end of NumTrials
            bestIndex = iga->get_BestOne();
            fprintf(ptr_EMAcalProfile,"%e %1.3f\n", PopFitness[bestIndex], obj->best_cc);
        }
    }
    for(i=0;i<ExecutionRuns;i++)
    {
        for(j=0;j<=NumOrigParam;j++)
        {
            if(j==NumOrigParam)
                printf("%e ",ModelParameters[i][j]);
            else
                printf("%f ",ModelParameters[i][j]);
        }
        printf("%1.3f\n",obj->best_cc);
    }
    
    if(ModelParameters)
    {
        for(i=0;i<ExecutionRuns;i++)
            delete [] ModelParameters[i];
        delete [] ModelParameters;
    }
    
    if(Chromosomes)
    {
        for(i=0;i<POP_SIZE;i++)
            delete [] Chromosomes[i];
        delete [] Chromosomes;
    }
    
    if(PopFitness)
        delete [] PopFitness;
    if(range)
        delete [] range;
    if(choose)
        delete [] choose;
    if(LocMask)
        delete [] LocMask;
    
    if(knowledge)
        delete knowledge;
    if(hill)
        delete hill;
    
    if(OA)
        delete OA;
    
    if(iga)
        delete iga;
    
    if(obj)
        delete obj;
     
    fclose(ptr_EMAcalProfile);
    fclose(ptr_Inputprofile);
    
    return 0;
}
void exit_with_help()
{
	printf(
           "Usage: HFODE_multiply [options] init_file kowledge_file\n"
           "options:\n"
           "-i : No.[?] of gene to run (Start from 0)\n"
           "-n : Number of run (default 30)\n"
           "-G : Number of generation (default 10000)\n"
           "-I : set the iteration of EMA (default 0)\n"
           "-F : Objective function (default 2)\n"
           "    0 -- (Xcal-Xexp)/Xexp\n"
           "    1 -- (Xcal-Xexp)\n"
           "    2 -- delta_diff + f0\n"
           "    3 -- delta_diff + f1\n"
           "    4 -- delta_diff\n"
           "-P : Number of popluation size (default 100)\n"
           "-t : Time points (necessary in simulated exp., default 10)\n"
           "-m change_mode : Test mode or RUN mode (default 1)\n"
           "	0 -- RUN mode\n"
           "	1 -- TEST mode\n"
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
	GeneIndex=0;
	ExecutionRuns=30;
	GENERATION=10000;
	NO_timepoints=10;
	POP_SIZE=100;
	FUNC=2;
	test_mode=1;
	cc_mode=0;
	iteration = 0;
	init_file_idx = 2;
	knowledge_idx = 3;
	i = 1;
        //======================//
    
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
			case 'n':
				ExecutionRuns = atoi(argv[i+1]);
				i++;
				break;
			case 'G':
				GENERATION = atoi(argv[i+1]);
				i++;
				break;
			case 'I':
				iteration = atoi(argv[i+1]);
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
        knowledge_idx = i;
    else{
        fprintf(stderr, "No knowledge file!!\n");
        exit_with_help();
    }
}
