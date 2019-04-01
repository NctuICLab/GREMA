#ifndef MODEL_H
#define MODEL_H
#include "Define.h"
#define REPLACE_ZERO		XMIN
//extern double AddGaussianNoise(double stdev, double mean);
struct Values
{
    double **X;
    double **dX;
    //double **fitness;
    //      double MError;
};
class Model
{
	//friend class S_System;
	friend struct Values;
public:
    Model();
	//S_System_Divide(int nVar, int nRun, int nTrial, double dt, bool doPerturb = false, int nPerturb = 0, double noise_deg = 0);
	~Model();
	int NumVariables() const { return NumVar; }
    int NumTFsVariables() const { return NumTF;}
	int NumRuns() const { return NumRun; }
	int NumTrials() const { return NumTrial; }
    int NumNochange() const { return Nochange;}
	bool DoPerturbation() const { return DoPerturb; }
	int NumPerturbationSet() const { return NumPerturb; }
    bool IsAddOriginalNoise() const { return IsAddOriNoise; }
    double OriginalNoiseDegree() const { return OriNoiseDegree; }
    bool IsSimulated() const { return isSimulated; }
    
    double ExpValue(int run, int i, int t) ;

    double CalValue(int run, int i, int t) ;

    //Set Hill function's parameters
	double BufValue(int run, int i, int t);
    void SetExpValues(int r, int i, int t, double val);
    void SetExpb(int i, double val);
    void SetExpTransMax(int i, double val);
    //void SetExpbeta(int i, int j, double val);
    void SetExpn(int i, int j, double val);
    void SetExpn(int i,int j, int val);
    void SetExpk(int i, int j, double val);
    void SetExpdegrade(int i, double val);
    
    void SetCalb(int i, double val);
    void SetCalTransMax(int i, double val);
    //void SetCalbeta(int i, int j, double val);
    //void SetCaln(int i, int j, int val);
    void SetCaln(int i, int j, double val);
    void SetCalk(int i, int j, double val);
    void SetCaldegrade(int i, double val);
    
	//set initial values
	int ReadInitFromFile(char* filename);
	void SetInitExpAllRun(double** X); //X[r][i][0]: run r, var i, trial 0
	void SetInitExpRun(int nRun, double *X);
	void SetInitCalAllRun(double** X); //X[r][i][0]: run r, var i, trial 0
	void SetInitCalRun(int nRun, double *X);
    
	//calcualte values
	//Experiment
	void CalAllExpValues();
	void CalRunExpValues(int nRun);
    
	//Calculated
	void CalAllCalValues(int index);
    void CalRunVarCalValues_Divide(int nRun, int index);//, int Var);
    
	//print
    void PrintExpModels(FILE* fout);
	void PrintAllExpXs(FILE* fout);
	void PrintExpXs(FILE* fout, int run);
	void PrintExpXs(FILE* fout, int run, int var);
	void PrintExpX(FILE* fout, int run, int var, int trial);
	void PrintCalModels(FILE* fout);
	void PrintRawCalModels(FILE* fout);
	void PrintAllCalXs(FILE* fout);
	void PrintCalXs(FILE* fout, int run);
	void PrintCalXs(FILE* fout, int run, int var);
	void PrintCalX(FILE* fout, int run, int var, int trial);
    
	void NormalizeExpValues();
    double AddGaussianNoise(double stdev, double mean);
	
    void SetDeltaS(double ds) { delta_s = ds; }
	double GetDeltaS() { return delta_s; }
    
	//void SetDeltaT(double dt) { delta_t = dt; }
	double GetDeltaT() { return delta_t; }
    //double GetsimDeltaT() { return sim_deltat;}
	double CalDeltaValue(int run, int i, int t);
	//double ExpDeltaValue(int run, int i, int t);
	//double CalFitnessValue(int run, int i, int t);
    void SetGeneIndex(int index){
        m_GeneIndex = index;
    }
    int b_upvalue() const { return b_up;}
    int b_downvalue() const { return b_down;}
    int TransMax_upvalue() const {return TransMax_up;}
    int TransMax_downvalue() const {return TransMax_down;}
    //int beta_upvalue() const { return beta_up;}
    //int beta_downvalue() const { return beta_down; }
    int n_upvalue() const { return n_up;}
    int n_downvalue() const { return n_down;}
    int k_upvalue() const { return k_up;}
    int k_downvalue() const { return k_down;}
    int degrade_upvalue() const { return degrade_up; }
    int degrade_downvalue() const { return degrade_down; }
    //double deltat_upvalue() { return deltat_up;}
    //double deltat_downvalue() { return deltat_down;}
    //int GetIndex(){return m_GeneIndex;}
private:
	void InitializeExp();
	void InitializeCal();
	void InitializeBufValues();
    //double deltat;
    int m_GeneIndex;
    //unsigned int CheckIndex;
	int NumVar; //num of total gene(TF+gene)
    int NumTF;// num of TFs
	int NumRun; //num of Runs
	int NumTrial; //num of trial(time point)
    int Nochange;
    //range of variables
    //double deltat_up;
    //double deltat_down;
    int b_up;
    int b_down;
    int TransMax_up;
    int TransMax_down;
    //int beta_up;
    //int beta_down;
    int n_up;
    int n_down;
    int k_up;
    int k_down;
    int degrade_up;
    int degrade_down;
    
	//for add perturbation set
	bool DoPerturb; //do perturbation or not
	int NumPerturb; //number of Perturbation set
	double NoiseDegree;
    
	//for add noise to observed data
	bool IsAddOriNoise;
	double OriNoiseDegree;
    
	bool* isSetInitExpVal;
	bool isInitExpModel;
    
	bool* isSetInitCalVal;
	bool isInitCalModel;
	
	double delta_t; //delta time
	double delta_s; //for structure skeletelizing
    
	bool isSimulated;
    //experiment value
    double* expb;
    double* expTransMax;
    //double**    exp_beta;
    double** exp_n;
    //int**       exp_n;
    double** exp_k;
    double* expdegrade;
    Values* expValues;
    //cal val
    double* calb;
    double* calTransMax;
    //double**    cal_beta;
    double** cal_n;
    //int**       cal_n;
    double** cal_k;
    double* caldegrade;
    Values* calValues;
    int count;

};
#endif 


