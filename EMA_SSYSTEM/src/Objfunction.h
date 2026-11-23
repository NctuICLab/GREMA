#ifndef Objfunction_h
#define Objfunction_h
#include "S-system-div.h"
#include "getknowledge.h"
#include "stdafx.h"
class Objfunction
{
public:
    Objfunction(S_System_Divide *ssystem,Knowledge *know);
    ~Objfunction();
    S_System_Divide *s_system;
    Knowledge *knowledge;
    void init_obj(int time_step,int type,int run);
    void Decode(double *genes);
    double fitness(double *genes);
    double fitnessValue(double *genes);
    int get_time_interval(){return time_interval;}
    //double get_EVAL(int i){return EVAL[i];}
    double  *EVAL;
    int	NumVars;
    int NumRuns;
    int NumTrials;
    int NumTFs;
    int PARAM_NUM;
    int NumGeneParam;
    int NO_POP;
    double  last_cc;//correlation coefficient of the last evaluated model (over all replicates)
    double  *RunCC;//[execution run] cc of the best model of every run, see RecordRun
    unsigned int    function_cost;
    double *GeneParamValues;
    double **ModelParameters;
    double  **cal_profile;
    void printCurrentModelParameters(double *genes,int type);
    //Re-evaluate the best chromosome of an execution run and keep its cc (RunCC[run]) and the
    //simulated profile of every replicate (cal_profile). function_cost is not changed.
    void RecordRun(int run,double *genes);
private:
    int g_NumConnections;
    int h_NumConnections;
    int time_interval;
    int runtime;
    int fitness_type;
};

#endif
