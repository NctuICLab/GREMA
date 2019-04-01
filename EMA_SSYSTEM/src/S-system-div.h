#ifndef S_SYSTEM_DIVIDE_H
#define S_SYSTEM_DIVIDE_H
#include "stdafx.h"
#define REPLACE_ZERO		XMIN
//extern double AddGaussianNoise(double stdev, double mean);
struct Values
{
    double **X;
    double **dX;
    double **fitness;
    //      double MError;
};
class S_System_Divide
{
	//friend class S_System;
	friend struct Values;
public:
	S_System_Divide();
	//S_System_Divide(int nVar, int nRun, int nTrial, double dt, bool doPerturb = false, int nPerturb = 0, double noise_deg = 0);
	~S_System_Divide();
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
    
	//experimental parameters and data in S-system
    /*
    double* ExpAlpha()   { return expAlpha; }
    double ExpAlpha(int i) { return expAlpha[i]; }
    double* ExpBeta()    { return expBeta; }
    double ExpBeta(int i)  { return expBeta[i]; }
    double** ExpH()	{ return exp_h; }
    double** ExpG()	{ return exp_g; }
    double ExpH(int i, int j) { return exp_h[i][j]; }
    double ExpG(int i, int j) { return exp_g[i][j];}
     */
    double ExpValue(int run, int i, int t) ;

	//calculated parameters and data in S-system
    /*
    double* CalAlpha()   { return calAlpha; }
    double CalAlpha(int i) { return calAlpha[i]; }
    double* CalBeta()    { return calBeta; }
    double CalBeta(int i)  { return calBeta[i]; }
    double** CalH()	{ return cal_h; }
    double** CalG()	{ return cal_g; }
    double CalH(int i, int j) { return cal_h[i][j]; }
    double CalG(int i, int j) { return cal_g[i][j];}
     */
    double CalValue(int run, int i, int t) ;

    //Set S-system's parameters
	double BufValue(int run, int i, int t);
    void SetExpValues(int r, int i, int t, double val);
    void SetExpAlpha(int i, double val);
    void SetExpBeta(int i, double val);
    void SetExpH(int i, int j, double val);
    void SetExpG(int i, int j, double val);
    
    void SetCalAlpha(int i, double val);
    void SetCalBeta(int i, double val);
    void SetCalH(int i, int j, double val);
    void SetCalG(int i, int j, double val);
    
	//set initial values
	int ReadInitFromFile(char* filename);
	void SetInitExpAllRun(double** X); //X[r][i][0]: run r, var i, trial 0
	void SetInitExpRun(int nRun, double *X);
	void SetInitCalAllRun(double** X); //X[r][i][0]: run r, var i, trial 0
	void SetInitCalRun(int nRun, double *X);
    
	//calcualte values
	//Experiment
	int CalAllExpValues();
	int CalRunExpValues(int nRun);
    
	//Calculated
	int CalAllCalValues();
    int CalRunVarCalValues_Divide(int nRun);//, int Var);
    
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
    
    
	
    
	//Buffer Value operations
	void CopyExpToBufValues();
	void CopyCalToBufValues();
	void SetBufValues(int run, double **X);
	void SetBufValues(int run, int var, double* X);
	void SetBufValues(int run, int var, int trial, double x);
    int CalAllCalValues_Divide_use_Buf();
	int CalRunCalValues_Divide_use_Buf(int nRun);
	int CalRunVarCalValues_Divide_use_Buf(int nRun, int Var);
    
	void NormalizeExpValues();
    double AddGaussianNoise(double stdev, double mean);
	
    void SetDeltaS(double ds) { delta_s = ds; }
	double GetDeltaS() { return delta_s; }
    
	void SetDeltaT(double dt) { delta_t = dt; }
	double GetDeltaT() { return delta_t; }
    double GetsimDeltaT() { return sim_deltat;}
	double CalDeltaValue(int run, int i, int t);
	double ExpDeltaValue(int run, int i, int t);
	double CalFitnessValue(int run, int i, int t);
	void SetGeneIndex(int index);
    double alpha_upvalue() { return alpha_up;}
    double alpha_downvalue() { return alpha_down;}
    double beta_upvalue() { return beta_up;}
    double beta_downvalue() { return beta_down;}
    double g_upvalue() { return g_up;}
    double g_downvalue() { return g_down;}
    double h_upvalue() { return h_up;}
    double h_downvalue() { return h_down;}
    double deltat_upvalue() { return deltat_up;}
    double deltat_downvalue() { return deltat_down;}
    int			m_GeneIndex;
private:
	void InitializeExp();
	void InitializeCal();
	void InitializeBufValues();
    double sim_deltat;
	int NumVar; //num of total gene(TF+gene)
    int NumTF;// num of TFs
	int NumRun; //num of Runs
	int NumTrial; //num of trial(time point)
    int Nochange;
    //range of variables
    double alpha_up;
    double alpha_down;
    double beta_up;
    double beta_down;
    double g_up;
    double g_down;
    double h_up;
    double h_down;
    double deltat_up;
    double deltat_down;
    
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
    double* expAlpha;
    double* expBeta;
    double** exp_g;
    double** exp_h;
    Values* expValues;
    //cal val
    double*		calAlpha;
    double*		calBeta;
    double**	cal_g;
    double**	cal_h;
    
    Values*		calValues;
    Values*		bufValues;

};
#endif //S_SYSTEM_DIVIDE_H


