#ifndef IGA_h
#define IGA_h

// Default values of IGA
#define DEFAULT_PC					0.8
#define DEFAULT_PM1					0.02
#define DEFAULT_PM2					0.05
#define DEFAULT_PS					0.2

#include "stdafx.h"
#include "S-system-div.h"
#include "getknowledge.h"
#include "Objfunction.h"
#include "OrthogonalArray.h"

class IGA
{
public:
    IGA(COrthogonalArray *OAtable,S_System_Divide *ssystem,Knowledge *know,Objfunction *objection);
    ~IGA();
    COrthogonalArray *OA;
    S_System_Divide *s_system;
    Knowledge *knowledge;
    Objfunction *obj;
    void    run(int run);
    int get_PARAM_NUM(){return PARAM_NUM;};
    int get_NumGeneParam(){return NumGeneParam;};
    int get_BestOne(){return BestOne;}
    void set_parameter(int GeneIndex,int iteration,int run);
    void initialIGA();
    void RangeVariable();
    void ChangeMutationRange();
    double MutateOneValue(double oldValue,int index,int fixloc);
    int NO_GEN;
    int NO_POP;
    int NO_g;
    int NO_h;
    int	NumVars, NumRuns, NumTrials, NumTFs,Num_nochange;
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
    void    initialpop();
    void    Evaluate();
    void    Selection();
    void    Crossover();
    void    Mutation();
    void    CX_mutation(unsigned int index_one,unsigned int index_two);
    unsigned int	*target;
    int 	NumCrossover;
    double	*Chromosomes;
    double  pc;
    double  ps;
    double  pm1;
    double  pm2;
    int m_size1,m_size2,s_size;
    int PARAM_NUM;
    int NumGeneParam;
    int *SortIndex;
    int alpha_index,beta_index,h_start;
    int BestOne;
    int gloc_possible;
    int hloc_possible;
    double  *Lower_bound,*Upper_bound,*Interval;
    bool same_interval;
    double *IntervalChangeByTemplate;
    double temperature,SA;
    int **Loc_tmp, **eachpop_loc, **sortloc;
    double  *BestValueEachGeneration;
    unsigned int    *NumFitnessFunction;
    unsigned int    cx1,cx2;
    int runtime;
    int runnow;
    int GeneIndex;
    //double  **cal_profile;
    FILE *ptr1;
    //---------OA----------------//
    int *DoMutationIndex;
    
};
#endif

