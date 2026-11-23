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
    //Mutation of a Gij/Hij/alpha/beta/delta_t value (fixloc: the gene of the Gij/Hij, negative = switched off).
    //Location slots are NOT mutated here, see MutateLocation.
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
    //---- gene slots -----------------------------------------------------------------------
    //w = 0: G (location,Gij) pairs at chromosome index 0.., w = 1: H pairs at h_start..
    //A location code is the gene index (gene is on) or -gene (neg_zero for gene 0, gene is off).
    //Slots 0..confirmed-1 are dedicated to the '+'/'-' genes and never change, so those
    //regulators are in every individual (also after crossover/mutation) with the sign of the knowledge.
    //The remaining slots hold the free '?' genes.
    int     SlotCount(int w){return (w==0)?NO_g:NO_h;}
    int     SlotBase(int w){return (w==0)?0:h_start;}
    int     ConfirmedCount(int w){return (w==0)?knowledge->get_gConfirm():knowledge->get_hConfirm();}
    int     PossibleCount(int w){return (w==0)?gloc_possible:hloc_possible;}
    bool    SlotsFixed(int w){return (w==0)?knowledge->get_fixg():knowledge->get_fixh();}
    bool    SlotsIdTF(int w){return (w==0)?knowledge->get_idgTF():knowledge->get_idhTF();}
    static int  LocGene(double code);
    static bool LocOn(double code){return code >= 0;}
    static double LocCode(int gene,bool on);
    void    BuildSortLoc();
    void    InitLocations(int w,double *chrom);
    double  RandomValue(int w,int gene,int index);
    double  MutateLocation(int w,const double *chrom,int slot);
    double  BestCC();
    void    RepairLocations(double *chrom);
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
    int **sortloc;
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

