#ifndef Objfunction_h
#define Objfunction_h
#include "Model.h"
#include "Getknowledge.h"
#include "Define.h"
class Objfunction
{
public:
    //Objfunction(Model *hill_f);
    Objfunction();
    ~Objfunction();
    void init_obj(int time_step,int type,int run, int cc, int No_reg, int index, Model *hill_f, double **Chrom);
    void Decode(int PopIndex);
    double fitness(int PopIndex);
    double fitnessValue(double *OneChrom);
    double OAfitenss(double **OATable, int OATableIndex);
    int get_time_interval(){return time_interval;}
    //double* get_EVAL(int i){
    //    double* fitness = &EVAL[i];
    //    return fitness;
    //}
    //void set_EVAL(int i, double fitnessValue){EVAL[i] = fitnessValue;}
    //double  *EVAL;
    int	NumVars;
    int NumRuns;
    int NumTrials;
    int NumTFs;
    int PARAM_NUM;
    int NumGeneParam;
    int NO_POP;
    double  best_cc;
    unsigned int    function_cost;
    double *GeneParamValues;
    double **Chromosomes;
    //double **ModelParameters;
    double  **cal_profile;
    void printCurrentModelParameters(int PopIndex);
private:
    Model *hill;
    int NumConnections;
    //int h_NumConnections;
    int time_interval;
    double min_fitness;
    int fitness_type;
    int runtime;
    int GeneIndex;
    int cc_mode;
    //double  *EVAL;
};
#endif
