#ifndef IGA_h
#define IGA_h

// Default values of IGA
/*
#define DEFAULT_PC					0.8
#define DEFAULT_PM1					0.02
#define DEFAULT_PM2					0.05
#define DEFAULT_PS					0.2
*/
#include "Define.h"
#include "Model.h"
#include "Getknowledge.h"
#include "Objfunction.h"
#include "OrthogonalArray.h"

class IGA
{
public:
    //IGA(COrthogonalArray *OAtable,Model *hill_f,Getknowledge *know,Objfunction *objection);
    IGA();
    ~IGA();
    
    void    run(int run);
    int get_PARAM_NUM(){return PARAM_NUM;};
    int get_NumGeneParam(){return NumGeneParam;};
    int get_BestOne(){return BestOne;}
    void set_parameter(int Index,int iteration, int run_mode,COrthogonalArray *OAtable,Model *hill_f,Getknowledge *know,Objfunction *objection,double **OrigChromosomes, double *fitnessvalue, int *Lock, int *Choose, int *Range, double **ModelP);
    void initialIGA();
    void initialpop();
    void RangeVariable();
    void ChangeMutationRange();
    double MutateOneValue(double oldValue,int index,int fixloc);
    
    int NO_GEN;
    int NO_POP;
    int NO_REG;
    int	NumVars, NumRuns, NumTrials, NumTFs,Num_nochange;
    int PARAM_NUM;
    int NumGeneParam;
    //double **ModelParameters;

    // return a random value in the interval [0, 1]
    inline double RandomFloat()
    {
        return (double)rand()/(RAND_MAX+1.0);
    }
    // return a random integer from 0 to n-1
    inline int RandomInt(int n)
    {
        return (int)(RandomFloat()*n);
    }
private:
    COrthogonalArray *OA;
    Model *hill;
    Getknowledge *knowledge;
    Objfunction *obj;
    void    Evaluate();
    void    Selection();
    void    Crossover();
    void    Mutation();
    void    CX_mutation(unsigned int index_one,unsigned int index_two);
    unsigned int	*target;
    int 	NumCrossover;
    double	**Chromosomes;
    //double **testChrom;
    double  pc;
    double  ps;
    double  pm1;
    double  pm2;
    int m_size1,m_size2,s_size;
    
    int loc_possible;
    int *SortIndex;
    int k_start,b_index,TransMax_index,degrade_index;
    int BestOne;
    double *PopFitness;
    double  *Lower_bound,*Upper_bound,*Interval;
    //bool same_interval;
    double *IntervalChangeByTemplate;
    double temperature,SA;
    int *LockMaskValue, *ChooseValue, *RangeValue;
    int *Loc_tmp, *eachpop_loc, *sortloc;
    double **ModelParameters;
    //double  *BestValueEachGeneration;
    //unsigned int    *NumFitnessFunction;
    unsigned int    cx1,cx2;
    int test_mode;
    //int runtime;
    int runnow;
    int Index;
    int MaxMask;
    //double  **cal_profile;
    FILE *ptr_progress;
    //---------OA----------------//
    int *DoMutationIndex;
    
};
#endif

