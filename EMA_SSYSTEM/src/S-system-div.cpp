
#include "S-system-div.h"
extern int ExecutionRuns;

S_System_Divide::S_System_Divide()//initialization
{
	NumVar = 0;
	NumRun = 0;
	NumTrial = 0;
    
	delta_t = 0;
	delta_s = 0;
    
	DoPerturb = false;
	NumPerturb = 0;
	NoiseDegree = 0;
    
	IsAddOriNoise = false;
	OriNoiseDegree = 0;
    
	isSimulated = false;
    
	//experiment
    expAlpha = NULL;
    expBeta = NULL;
    exp_g = NULL;
    exp_h = NULL;
	
	isSetInitExpVal = NULL;
	expValues = NULL;
	isInitExpModel = false;
    
	//calculated
    calAlpha = NULL;
    calBeta = NULL;
    cal_g = NULL;
    cal_h = NULL;
    
    isSetInitCalVal = NULL;
	calValues = NULL;
	isInitCalModel = false;
	bufValues = NULL;
	m_GeneIndex=0;
    //printf("finish init s\n");
}
/*
S_System_Divide::S_System_Divide(int nVar, int nRun, int nTrial,double dt, bool doPerturb, int nPerturb, double noise_deg)
: NumVar(nVar), NumRun(nRun), NumTrial(nTrial), delta_t(dt), DoPerturb(doPerturb), NumPerturb(nPerturb), NoiseDegree(noise_deg)
{
	InitializeExp();
	InitializeCal();
	bufValues = NULL;
}
*/
S_System_Divide::~S_System_Divide()
{
	int i, j;
    if(NumVar != 0)
    {
        //experiment
        if(expAlpha) delete [] expAlpha;
        if(expBeta)  delete [] expBeta;
        if(exp_g)
        {
            for(i = 0; i < NumVar; i++)
                delete [] exp_g[i];
            delete [] exp_g;
        }
        
        if(exp_h)
        {
            for(i = 0; i < NumVar; i++)
                delete [] exp_h[i];
            delete [] exp_h;
        }
        //calculated
        if(calAlpha) delete [] calAlpha;
        if(calBeta)  delete [] calBeta;
        if(cal_g)
        {
            for(i = 0; i < NumVar; i++)
                delete [] cal_g[i];
            delete [] cal_g;
        }
        if(cal_h)
        {
            for(i = 0; i < NumVar; i++)
                delete [] cal_h[i];
            delete [] cal_h;
        }
    }
  
	//experiment
	if(expValues != NULL)
	{
		for(i = 0; i < NumRun; i++)
		{
			for(j = 0; j < NumVar; j++)
			{
                delete [] expValues[i].X[j];
				delete [] expValues[i].dX[j];
				delete [] expValues[i].fitness[j];
			}
            delete [] expValues[i].X;
			delete [] expValues[i].dX;
			delete [] expValues[i].fitness;
		}
		delete [] expValues;
	}
    
	//calculated
	if(calValues != NULL)
	{
		for(i = 0; i < NumRun; i++)
		{
			for(j = 0; j < NumVar; j++)
			{
				delete [] calValues[i].X[j];
				delete [] calValues[i].dX[j];
                delete [] calValues[i].fitness[j];
			}
			delete [] calValues[i].X;
			delete [] calValues[i].dX;
            delete [] calValues[i].fitness[j];
		}
		delete [] calValues;
	}

	//Experiment
	if(isSetInitExpVal != NULL)
		delete [] isSetInitExpVal;
    
	//calculated
	if(isSetInitCalVal != NULL)
		delete [] isSetInitCalVal;
    
	//Buffer
	if(bufValues != NULL)
	{
		for(i = 0; i < NumRun; i++)
		{
			for(j = 0; j < NumVar; j++)
			{
				//delete [] bufValues[i].exp[j];
				//delete [] bufValues[i].cal[j];
				//delete [] bufValues[i].delta[j];
				
				delete [] bufValues[i].X[j];
				delete [] bufValues[i].dX[j];
			}
			//delete [] bufValues[i].exp;
			//delete [] bufValues[i].cal;
			//delete [] bufValues[i].delta;
			
			delete [] bufValues[i].X;
			delete [] bufValues[i].dX;
		}
		delete [] bufValues;
	}
}

void S_System_Divide::InitializeExp()
{
	int i, j;
    //Allocate memory space for S-system model
    expAlpha = new double[NumVar];
    expBeta = new double[NumVar];
    
    memset(expAlpha, 0, NumVar*sizeof(double));
    memset(expBeta, 0, NumVar*sizeof(double));
    
    exp_g = new double* [NumVar];
    exp_h = new double* [NumVar];
    for(i = 0; i < NumVar; i++)
    {
        exp_g[i] = new double[NumVar];
        exp_h[i] = new double[NumVar];
        memset(exp_g[i], 0, NumVar*sizeof(double));
        memset(exp_h[i], 0, NumVar*sizeof(double));
    }
    
	//Allocate memory space for Values
	expValues = new Values[NumRun+NumPerturb];
	for(i = 0; i < NumRun+NumPerturb; i++)
	{
		expValues[i].X = new double* [NumVar];
		expValues[i].dX = new double* [NumVar];
		expValues[i].fitness = new double* [NumVar];
        
		for(j = 0; j < NumVar; j++)
		{
			expValues[i].X[j] = new double[NumTrial+1];
			expValues[i].dX[j] = new double[NumTrial+1];
			expValues[i].fitness[j] = new double[NumTrial+1];
		}
	}
    //printf("InitializeExp3 NumVar=%d NumTF=%d NumRun=%d\n",NumVar,NumTF,NumRun);
	isSetInitExpVal = new bool[NumRun+NumPerturb];
	for(i = 0; i < NumRun+NumPerturb; i++)
		isSetInitExpVal[i] = false;
	
	isInitExpModel = true;
    //printf("InitializeExp over NumVar=%d NumTF=%d NumRun=%d\n",NumVar,NumTF,NumRun);
}

void S_System_Divide::InitializeCal()
{
	int i, j;
	//Allocate memory space for new model
    //Allocate memory space for S-system model
    calAlpha = new double[NumVar];
    calBeta = new double[NumVar];
    
    cal_g = new double* [NumVar];
    cal_h = new double* [NumVar];
    for(i = 0; i < NumVar; i++)
    {
        cal_g[i] = new double[NumVar];
        cal_h[i] = new double[NumVar];
    }
	//Allocate memory space for Values
	calValues = new Values[NumRun+NumPerturb];
	for(i = 0; i < NumRun+NumPerturb; i++)
	{
		calValues[i].X = new double* [NumVar];
		calValues[i].dX = new double* [NumVar];
		calValues[i].fitness = new double* [NumVar];
        
		for(j = 0; j < NumVar; j++)
		{
			calValues[i].X[j] = new double[NumTrial+1];
			calValues[i].dX[j] = new double[NumTrial+1];
			calValues[i].fitness[j] = new double[NumTrial+1];
		}
	}
    
    
	isSetInitCalVal = new bool[NumRun+NumPerturb];
	for(i = 0; i < NumRun+NumPerturb; i++)
		isSetInitCalVal[i] = false;
    
	isInitCalModel = true;
}
double S_System_Divide::AddGaussianNoise(double stdev, double mean)
{
    double r1, r2;
    double x;
    
    do
    {
        r1 = rand() / (double)RAND_MAX;
    }while(r1 == 0);
    
    r2 = rand() / (double)RAND_MAX;
    
    x = (stdev * sqrt( -2.0*log(r1) ) * cos(2*PI*r2) ) + mean;
    
    return x;
}

int S_System_Divide::ReadInitFromFile(char* filename)
{
	int rtn = 1;
	int i, j, r, run, t;
	FILE* fp;
	char line[BUF_SIZE];
	char* token;
	double** X;
	double tmp;
	fp = fopen(filename, "r");
	if(fp == NULL)
	{
		fprintf(stderr, "%s open failed!!!\n", filename);
		return -1;
	}
    
	while(fgets(line, BUF_SIZE-1, fp))
	{
		if(line[0] == '#') //skip comment line
			continue;//leave while loop
		
		//skip separate characters
		token = strtok(line, "\t\r\n #"); 
		if(token == NULL)
			continue;
        
		if(strlen(token) == 0)
			continue;
        
		if(strcmp(token, "num_var") == 0)
		{
			token = strtok(NULL, "\t\r\n #");
			NumVar = atoi(token);
            //printf("NumVar=%d\n",NumVar);
        }
        else if(strcmp(token, "num_tf") == 0)
		{
			token = strtok(NULL, "\t\r\n #");
			NumTF = atoi(token);
            //printf("NumTF=%d\n",NumTF);
		}
		else if(strcmp(token, "num_run") == 0)
		{
			token = strtok(NULL, "\t\r\n #");
			NumRun = atoi(token);
            //printf("NumRun=%d\n",NumRun);
		}
		else if(strcmp(token, "num_trial") == 0)
		{
			token = strtok(NULL, "\t\r\n #");
			NumTrial = atoi(token);
            //printf("NumTrial=%d\n",NumTrial);
		}
		else if(strcmp(token, "sim_deltat") == 0)
		{
			token = strtok(NULL, "\t\r\n #");
			sim_deltat = atof(token);
            //printf("sim_deltat=%f\n",sim_deltat);
		}
        
        else if(strcmp(token, "delta_s") == 0)
		{
			token = strtok(NULL, "\t\r\n #");
			delta_s = atof(token);
            //printf("delta_t=%f\n",delta_t);
		}
        else if(strcmp(token, "nochange") == 0)
        {
            token = strtok(NULL, "\t\r\n #");
            Nochange = atoi(token);
            //printf("Nochange=%d\n",Nochange);
        }
        else if(strcmp(token, "alpha_up") == 0)
        {
            token = strtok(NULL, "\t\r\n #");
            alpha_up = atoi(token);
            //printf("alpha_up=%f\n",alpha_up);
        }
        else if(strcmp(token, "alpha_down") == 0)
        {
            token = strtok(NULL, "\t\r\n #");
            alpha_down = atoi(token);
            //printf("alpha_down=%f\n",alpha_down);
        }
        else if(strcmp(token, "beta_up") == 0)
        {
            token = strtok(NULL, "\t\r\n #");
            beta_up = atoi(token);
            //printf("beta_up=%f\n",beta_up);
        }
        else if(strcmp(token, "beta_down") == 0)
        {
            token = strtok(NULL, "\t\r\n #");
            beta_down = atoi(token);
            //printf("beta_down=%f\n",beta_down);
        }
        else if(strcmp(token, "g_up") == 0)
        {
            token = strtok(NULL, "\t\r\n #");
            g_up = atoi(token);
            //printf("g_up=%f\n",g_up);
        }
        else if(strcmp(token, "g_down") == 0)
        {
            token = strtok(NULL, "\t\r\n #");
            g_down = atoi(token);
            //printf("g_downn=%f\n",g_down);
        }
        else if(strcmp(token, "h_up") == 0)
        {
            token = strtok(NULL, "\t\r\n #");
            h_up = atoi(token);
            //printf("h_up=%f\n",h_up);
        }
        else if(strcmp(token, "h_down") == 0)
        {
            token = strtok(NULL, "\t\r\n #");
            h_down = atoi(token);
            //printf("h_downn=%f\n",h_down);
        }
        else if(strcmp(token, "deltat_up") == 0)
        {
            token = strtok(NULL, "\t\r\n #");
            deltat_up = atof(token);
            //printf("deltat_up=%f\n",deltat_up);
        }
        else if(strcmp(token, "deltat_down") == 0)
        {
            token = strtok(NULL, "\t\r\n #");
            deltat_down = atof(token);
            //printf("deltat_down=%f\n",deltat_down);
        }
		else if(strcmp(token, "add_perturbation") == 0)
		{
			token = strtok(NULL, "\t\r\n #");
			//for(i = 0; i < (int)strlen(token); i++)
			//	token[i] = tolower(token[i]);
            
			if(strcmp(token, "yes") == 0)
			{
				DoPerturb = true;
				
				token = strtok(NULL, "\t\r\n #");
				NumPerturb = atoi(token);
                
				token = strtok(NULL, "\t\r\n #");
				NoiseDegree = atof(token);
			}
			else
			{
				DoPerturb = false;
				NumPerturb = 0;
				NoiseDegree = 0;
			}
		}
		else if(strcmp(token, "ori_noise") == 0)
		{
			token = strtok(NULL, "\t\r\n #");
			if(strcmp(token, "yes") == 0)
			{
				IsAddOriNoise = true;
                
				token = strtok(NULL, "\t\r\n #");
				OriNoiseDegree = atof(token);
			}
			else
			{
				IsAddOriNoise = false;
				OriNoiseDegree = 0;
			}
		}
		else if(strcmp(token, "init_data") == 0)//initial data(1st time point data)
		{
            //printf("init data\n");
			if( (NumVar > 0 && (NumRun+NumPerturb) > 0 && NumTrial > 0) )
			{
				if(!isInitExpModel)
					InitializeExp();
				if(!isInitCalModel)
					InitializeCal();
			}
			else
			{
				fprintf(stderr, "Error input data file, initalization of data failed!!!\n");
				exit(-1);
			}
            
			//allocate memory space
			X = new double* [NumRun];
			for(i = 0; i < NumRun; i++)
				X[i] = new double[NumVar];
            
			//Read init. data
			for(i = 0; i < NumRun; i++)
			{
				fgets(line, BUF_SIZE - 1, fp);
                token = strtok(line, " \t\r\n");
				for(j = 0; j < NumVar; j++)
				{
					X[i][j] = atof(token);
					token = strtok(NULL, " \t\r\n");
				}
			}
            
			SetInitExpAllRun(X);
			SetInitCalAllRun(X);
            
			for(i = 0; i < NumRun; i++)
				delete [] X[i];
			delete [] X;
		}
		else if(strcmp(token, "init_parameters") == 0)
		{
            //printf("init parameters\n");
			if( (NumVar > 0 && NumRun > 0 && NumTrial > 0) )
			{
				if(!isInitExpModel)
					InitializeExp();
				if(!isInitCalModel)
					InitializeCal();
			}
			else
			{
				fprintf(stderr, "Error input data file, initalization of parameters failed!!!\n");
				exit(-1);
			}
            
			for(i = 0; i < NumVar; i++)
			{
				fgets(line, BUF_SIZE-1, fp);
                //alpha
                token = strtok(line, "\n\r\t #");
                expAlpha[i] = atof(token);
                //printf("expAlpha[%d]=%f\n",i,expAlpha[i]);
                //gij
                for(j = 0; j < NumVar; j++)
                {
                    token = strtok(NULL, "\n\r\t #");
                    exp_g[i][j] = atof(token);
                    //printf("exp_g[%d][%d]=%f\n",i,j,exp_g[i][j]);
                }
                
                //hij
                token = strtok(NULL, "\n\r\t #");
                expBeta[i] = atof(token);
                for(j = 0; j < NumVar; j++)
                {
                    token = strtok(NULL, "\n\r\t #");
                    exp_h[i][j] = atof(token);
                    //printf("exp_h[%d][%d]=%f\n",i,j,exp_h[i][j]);
                }

            }
			isSimulated = true;
		}
        
        
		else if(strcmp(token, "list_data") == 0)
		{
            
			if( (NumVar > 0 && NumRun > 0 && NumTrial > 0) )
			{
				if(!isInitExpModel)
					InitializeExp();
				if(!isInitCalModel)
					InitializeCal();
			}
			else
			{
				fprintf(stderr, "Error input data file, initalization of parameters failed!!!\n");
				exit(-1);
			}
			//for each run
			for(r = 0; r < NumRun; r++)
			{
				fgets(line, BUF_SIZE-1, fp);
				token = strtok(line, "\n\r\t #"); //id of r-th run
				run = atoi(token);
				//for each variable
				for(i = 0; i < NumVar; i++)
				{
					
                    fgets(line, BUF_SIZE-1, fp);
                    token = strtok(line, "\n\r\t #"); //exprssion value at time 0
                    tmp = atof(token);
					if(tmp == 0)
                    {
                        expValues[r].X[i][0] = XMIN;
                    }
                    else
                    {
                        expValues[r].X[i][0] = tmp;
                    }
					//calValues[r].X[i][0] = atof(token);
                    
					//for each trial (1 ~ NumTrial)
					for(j = 1; j < NumTrial; j++)
					{
						token = strtok(NULL, "\n\r\t #"); //exprssion value at time j
                        tmp = atof(token);
                        if(tmp == 0)
                        {
                            expValues[r].X[i][j] = XMIN;
                        }
                        else
                        {
                            expValues[r].X[i][j] = tmp;
                        }
					}
				}
                //printf("input gene expression\n");
				isSetInitExpVal[r] = true;
				//isSetInitCalVal[r] = true;
			}
            
			if(DoPerturb)
			{
				j = 0;
				for(r=NumRun; r < NumRun+NumPerturb; r++)// r: 幾組實驗
				{
					//perturbate original data to generate the data of perturbation
					for(t = 0; t < NumTrial; t++) //for each trial    每組取樣次數
					{
						for(i = 0; i < NumVar; i++) //for each variable
						{
							tmp = expValues[j].X[i][t];
							expValues[r].X[i][t] = AddGaussianNoise(tmp * NoiseDegree, tmp);
						}
					}
					
					j ++; //change original run for generating perturbation set
					j %= NumRun; //keep the idx in the range of [0, NumRun)
					
					isSetInitExpVal[r] = true;
				}
			}//end of DoPerturb
		}
		else
			continue;
	}
	fclose(fp);
    
	return rtn;
}
/*
 void S_System_Divide::PrintExpModels(FILE* fout)
 {
 int i, j;
 for(i = 0; i < NumVar; i++)
 {
 fprintf(fout, "%g", expAlpha[i]);
 for(j = 0; j < NumVar; j++)
 fprintf(fout, " %g", exp_g[i][j]);
 
 fprintf(fout, " %g", expBeta[i]);
 for(j = 0; j < NumVar; j++)
 fprintf(fout, " %g", exp_h[i][j]);
 fprintf(fout, "\n");
 }
 }
 */
void S_System_Divide::PrintAllExpXs(FILE* fout)
{
	int r;
	int cnt = NumRun + NumPerturb;
	for(r = 0; r < cnt; r++)
	{
		PrintExpXs(fout, r);
		fprintf(fout, "\n\n");
	}
}

void S_System_Divide::PrintExpXs(FILE* fout, int run)
{
	int cnt = NumRun + NumPerturb;
	if(run >= cnt)
	{
		fprintf(stderr, "Error in PrintExpXs: Wrong run range!!!\n");
		return;
	}
    
	int v;
	for(v = 0; v < NumVar; v++)
	{
		PrintExpXs(fout, run, v);
		fprintf(fout, "\n");
	}
}

void S_System_Divide::PrintExpXs(FILE* fout, int run, int var)
{
	int cnt = NumRun + NumPerturb;
	if(run >= cnt)
	{
		fprintf(stderr, "Error in PrintExpXs: Wrong run range!!!\n");
		return;
	}
    
	if(var >= NumVar)
	{
		fprintf(stderr, "Error in PrintExpXs: Wrong variable range!!!\n");
		return;
	}
    
	int t;
	for(t = 0; t < NumTrial; t++)
		PrintExpX(fout, run, var, t);
	//fprintf(fout, "\n");
}

void S_System_Divide::PrintExpX(FILE* fout, int run, int var, int trial)
{
	int cnt = NumRun + NumPerturb;
	if(run >= cnt)
	{
		fprintf(stderr, "Error in PrintExpX: Wrong run range!!!\n");
		return;
	}
    
	if(var >= NumVar)
	{
		fprintf(stderr, "Error in PrintExpX: Wrong variable range!!!\n");
		return;
	}
    
	if(trial >= NumTrial)
	{
		fprintf(stderr, "Error in PrintExpX: Wrong trial range!!!\n");
		return;
	}
    //printf("%f ",expValues[run].X[var][trial]);
	fprintf(fout, "%g\t", expValues[run].X[var][trial]);
}

void S_System_Divide::SetInitExpAllRun(double** X) //X[r][i][0]: run r, var i, trial 0
{
	int i;
	for(i = 0; i < NumRun; i++)
		SetInitExpRun(i, X[i]);
}

void S_System_Divide::SetInitExpRun(int nRun, double* X)
{
	int i;
	for(i = 0; i < NumVar; i++)
		expValues[nRun].X[i][0] = X[i];
	
	isSetInitExpVal[nRun] = true;
}

int S_System_Divide::CalAllExpValues()//simulated exp
{
	if(!isSimulated)
		return -1;
    
	int i, v, t, j;
	double tmp;
	int rtn = 1;
	for(i = 0; i < NumRun; i++)
	{
		rtn *= CalRunExpValues(i);
	}
    
	if(IsAddOriNoise)
	{
		for(i = 0; i < NumRun; i++)
		{
			//Add noise to the original observed data
			for(t = 0; t < NumTrial; t++) //for each trial
			{
				for(v = 0; v < NumVar; v++)//for each variable
				{
					tmp = expValues[i].X[v][t];
					expValues[i].X[v][t] = AddGaussianNoise(tmp * OriNoiseDegree, tmp);
				}//end of for each variable
			}//end of for each trial
		}
	}
    
	if(DoPerturb)
	{
		j = 0;
		for(i=NumRun; i < NumRun+NumPerturb; i++)
		{
			//perturbate original data to generate the data of perturbation
			for(t = 0; t < NumTrial; t++) //for each trial
			{
				for(v = 0; v < NumVar; v++) //for each variable
				{
					tmp = expValues[j].X[v][t];
                    expValues[i].X[v][t] = AddGaussianNoise(tmp * NoiseDegree, tmp);
				}
			}
			j ++; //change original run for generating perturbation set
			j %= NumRun; //keep the idx in the range of [0, NumRun)
		}
	}
    
	return rtn;
}

int S_System_Divide::CalRunExpValues(int nRun)//produce simulated profile
{
    double termA, termB, tmp, logval;
    int i, j, t;
    //double expval, val, X;
    
    if(isSetInitExpVal[nRun] != true)
    {
        fprintf(stderr, "Experiment Run %d is not initialized!\n", nRun);
        return -1;
    }
    
    for(t = 0; t < NumTrial; t++)
    {
        for(i = 0; i < NumVar; i++)
        {
            
            termA = 1;
            termB = 1;
            for(j = 0; j < NumVar; j++)
            {
                logval = log(expValues[nRun].X[j][t]);
                termA *= exp( exp_g[i][j] * logval );
                termB *= exp( exp_h[i][j] * logval );
            }
            
            termA *= expAlpha[i];
            termB *= expBeta[i];
            
            expValues[nRun].dX[i][t] = termA - termB;
            
            tmp = expValues[nRun].dX[i][t] * sim_deltat + expValues[nRun].X[i][t];
            if(tmp >= XMAX)
                expValues[nRun].X[i][t+1] = XMAX;
            else if(tmp <= XMIN)
                expValues[nRun].X[i][t+1] = XMIN;
            else//*/
                expValues[nRun].X[i][t+1] = tmp;
            
        }
    }//end of number of trial
    
    return 1;
}

void S_System_Divide::SetExpValues(int r, int i, int t, double val)
{
	if(r >= NumRun + NumPerturb)
	{
		fprintf(stderr, "Error: run(r), SetExpValues(r, i, t, val) out of range!!!\n");
		return;
	}
    
	if(i >= NumVar)
	{
		fprintf(stderr, "Error: variable(i), SetExpValues(r, i, t, val) out of range!!!\n");
		return;
	}
    
	if(t >= NumTrial)
	{
		fprintf(stderr, "Error: time point(t), SetExpValues(r, i, t, val) out of range!!!\n");
		return;
	}
    
	expValues[r].X[i][t] = val;
}

//--------------------------------------------------------------------------

void S_System_Divide::SetExpAlpha(int i, double val)
{
    if(i < NumVar && expAlpha != NULL)
        expAlpha[i] = val;
    else
        fprintf(stderr, "Error: SetExpAlpha out of index range!!!\n");
}

void S_System_Divide::SetExpBeta(int i, double val)
{
    if(i < NumVar && expBeta != NULL)
        expBeta[i] = val;
    else
        fprintf(stderr, "Error: SetExpBeta out of index range!!!\n");
}

void S_System_Divide::SetExpH(int i, int j, double val)
{
    if(i < NumVar && j < NumVar && exp_h != NULL)
        exp_h[i][j] = val;
    else
        fprintf(stderr, "Error: SetExpH out of index range!!!\n");
}

void S_System_Divide::SetExpG(int i, int j, double val)
{
    if(i < NumVar && j < NumVar && exp_g != NULL)
        exp_g[i][j] = val;
    else
        fprintf(stderr, "Error: SetExpG out of index range!!!\n");
}

void S_System_Divide::SetInitCalAllRun(double** X) //X[r][i][0]: run r, var i, trial 0
{
	int i;
	for(i = 0; i < NumRun; i++)
		SetInitCalRun(i, X[i]);
}

void S_System_Divide::SetInitCalRun(int nRun, double* X)
{
	int i;
	for(i = 0; i < NumVar; i++)
		calValues[nRun].X[i][0] = X[i];
	
	isSetInitCalVal[nRun] = true;
}

/*
void S_System_Divide::PrintCalModels(FILE* fout)
{
	int i, j;
	for(i = 0; i < NumVar; i++)
	{
		fprintf(fout, "%g", calAlpha[i]);
		for(j = 0; j < NumVar; j++)
			fprintf(fout, " %g", (fabs(cal_g[i][j]) <= delta_s ? 0 : cal_g[i][j]) );//if(fabs(cal_g[i][j]<=delta_s)){return 0;}else{return cal_g[i][j]} 
		fprintf(fout, " %g", calBeta[i]);
		for(j = 0; j < NumVar; j++)
			fprintf(fout, " %g", (fabs(cal_h[i][j]) <= delta_s ? 0 : cal_h[i][j]) );
		fprintf(fout, "\n");
	}
}

void S_System_Divide::PrintRawCalModels(FILE* fout)
{
	int i, j;
	for(i = 0; i < NumVar; i++)
	{
		fprintf(fout, "%g", calAlpha[i]);
		for(j = 0; j < NumVar; j++)
			fprintf(fout, " %g", cal_g[i][j] );
		
		fprintf(fout, " %g", calBeta[i]);
		for(j = 0; j < NumVar; j++)
			fprintf(fout, " %g", cal_h[i][j] );
		fprintf(fout, "\n");
	}
}
*/
void S_System_Divide::PrintAllCalXs(FILE* fout)
{
	int r;
	int cnt = ExecutionRuns;
	
	for(r = 0; r < cnt; r++)
	{
		PrintCalXs(fout, r);
        //fprintf(fout, "\n\n");
	}
}

void S_System_Divide::PrintCalXs(FILE* fout, int run)
{
	int cnt = ExecutionRuns;
	if(run >= cnt)
	{
		fprintf(stderr, "Error in PrintCalXs: Wrong run range!!!\n");
		return;
	}
    
	int v;
	for(v = 0; v < 1; v++)
	{
		PrintCalXs(fout, run, v);
		//fprintf(fout, "\n");
	}
}

void S_System_Divide::PrintCalXs(FILE* fout, int run, int var)
{
	int cnt = ExecutionRuns;
	if(run >= cnt)
	{
		fprintf(stderr, "Error in PrintCalXs: Wrong run range!!!\n");
		return;
	}
    
	if(var >= 1)
	{
		fprintf(stderr, "Error in PrintCalXs: Wrong variable range!!!\n");
		return;
	}
    
	int t;
	for(t = 0; t < NumTrial; t++)
		PrintCalX(fout, run, var, t);
	fprintf(fout, "\n");
}

void S_System_Divide::PrintCalX(FILE* fout, int run, int var, int trial)
{
	int cnt = ExecutionRuns;
	if(run >= cnt)
	{
		fprintf(stderr, "Error in PrintCalX: Wrong run range!!!\n");
		return;
	}
    
	if(var >= 1)
	{
		fprintf(stderr, "Error in PrintCalX: Wrong variable range!!!\n");
		return;
	}
    
	if(trial >= NumTrial)
	{
		fprintf(stderr, "Error in PrintCalX: Wrong trial range!!!\n");
		return;
	}
    
	//fprintf(fout, "%g ", calValues_temp[run].X[var][trial]);
}

int S_System_Divide::CalAllCalValues()//calculate the fitting profile
{
	int i, v;
	int rtn = 1;
    int totalnumRun = NumRun + NumPerturb;
	for(i = 0; i < totalnumRun; i++)	//i 組實驗
	{
		//if using real data, set the initial calculated value as the observed one
		if(!isSimulated)
		{
			if(isSetInitCalVal[i] != true)
			{
				//copy the perturbated result from expValues
				for(v = 0; v < NumVar; v++)
					calValues[i].X[v][0] = expValues[i].X[v][0];
                    				
				isSetInitCalVal[i] = true;
			}
		}
        rtn *=CalRunVarCalValues_Divide(i);//, j);	//run times, number of xi
	}
    return rtn;
}

void S_System_Divide::SetCalAlpha(int i, double val)
{
    if(i < NumVar && calAlpha != NULL)
        calAlpha[i] = val;
    else
        fprintf(stderr, "Error: SetCalAlpha out of index range!!!\n");
}

void S_System_Divide::SetCalBeta(int i, double val)
{
    if(i < NumVar && calBeta != NULL)
        calBeta[i] = val;
    else
        fprintf(stderr, "Error: SetCalBeta out of index range!!!\n");
}

void S_System_Divide::SetCalH(int i, int j, double val)
{
    if(i < NumVar && j < NumVar && cal_h != NULL)
        cal_h[i][j] = val;
    else
        fprintf(stderr, "Error: SetCalH out of index range!!!\n");
}

void S_System_Divide::SetCalG(int i, int j, double val)
{
    if(i < NumVar && j < NumVar && cal_g != NULL)
        cal_g[i][j] = val;
    else
        fprintf(stderr, "Error: SetCalG out of index range!!!\n");
}


int S_System_Divide::CalRunVarCalValues_Divide(int nRun)// calculate the fitting profile
{
    double termA, termB, tmp, logval;
    int i, j, t;
    int flag;
    int Var;
    //if using real data, set the initial calculated value as the observed one
    if(!isSimulated)
    {
        if(isSetInitCalVal[nRun] != true)
        {
            //copy the perturbated result from expValues
            for(j = 0; j < NumVar; j++)
                calValues[nRun].X[j][0] = expValues[nRun].X[j][0];
            
            isSetInitCalVal[nRun] = true;
        }
    }
    
    if(DoPerturb)
    {
        if(isSetInitCalVal[nRun] != true)
        {
            //copy the perturbated result from expValues
            for(i = 0; i < NumVar; i++)
                calValues[nRun].X[i][0] = expValues[nRun].X[i][0];
            
            isSetInitCalVal[nRun] = true;
        }
    }
    
    if(isSetInitCalVal[nRun] != true)
    {
        fprintf(stderr, "Calculated_Divide Run %d is not initialized!\n", nRun);
        return -1;
    }
    
    
    for(t = 0; t < NumTrial; t++)//第幾個測試
    {
        Var=m_GeneIndex;
        termA = 1;
        termB = 1;
        for(j = 0; j < NumVar; j++)	//	j第幾個參數
        {
            //Divide and conquer technique
            flag=0;
            if(j ==m_GeneIndex)
                flag=1;
            if(flag)
                logval = log(calValues[nRun].X[j][t]);
            else //if(i != j)
                logval = log(expValues[nRun].X[j][t]);
            termA *= exp(  cal_g[Var][j] * logval );
            termB *= exp(  cal_h[Var][j] * logval );
        }
            
        termA *= calAlpha[Var];
        termB *= calBeta[Var];
        calValues[nRun].dX[Var][t] = termA - termB;
        tmp = calValues[nRun].dX[Var][t] * delta_t + calValues[nRun].X[Var][t];
        //printf("%d tmp=%f\n",t,tmp);
        if(tmp >= XMAX)
            calValues[nRun].X[Var][t+1] = XMAX;
        else if(tmp <= XMIN)
            calValues[nRun].X[Var][t+1] = XMIN;
        else
            calValues[nRun].X[Var][t+1] = tmp;
    }//end of number of trial
    return 1;
}

//Buffer Value operations
void S_System_Divide::InitializeBufValues()
{
	int i, j;
	int cnt = NumRun + NumPerturb;
    
	//Allocate memory space for Values
	bufValues = new Values[cnt];
	for(i = 0; i < cnt; i++)
	{
		//bufValues[i].exp = new double* [NumVar];
		//bufValues[i].cal = new double* [NumVar];
		//bufValues[i].delta = new double* [NumVar];
		bufValues[i].X = new double* [NumVar];
		bufValues[i].dX = new double* [NumVar];
        
		for(j = 0; j < NumVar; j++)
		{
			//bufValues[i].exp[j] = new double[NumTrial];
			//bufValues[i].cal[j] = new double[NumTrial];
			//bufValues[i].delta[j] = new double[NumTrial];
			
			bufValues[i].X[j] = new double[NumTrial];
			bufValues[i].dX[j] = new double[NumTrial];
		}
	}
}

void S_System_Divide::CopyExpToBufValues()
{
	if(bufValues == NULL)
		InitializeBufValues();
    
	int i, j;
	int cnt = NumRun + NumPerturb;
	for(i = 0; i < cnt; i++)
	{
		for(j = 0; j < NumVar; j++)
		{
			memcpy(bufValues[i].X[j], expValues[i].X[j], sizeof(double)*NumTrial );
			memcpy(bufValues[i].dX[j], expValues[i].dX[j], sizeof(double)*NumTrial );
		}
	}
}

void S_System_Divide::CopyCalToBufValues()
{
	if(bufValues == NULL)
		InitializeBufValues();
    
	int i, j;
	int cnt = NumRun + NumPerturb;
	for(i = 0; i < cnt; i++)
	{
		for(j = 0; j < NumVar; j++)
		{
			memcpy(bufValues[i].X[j], calValues[i].X[j], sizeof(double)*NumTrial );
			memcpy(bufValues[i].dX[j], calValues[i].dX[j], sizeof(double)*NumTrial );
		}
	}
}

void S_System_Divide::SetBufValues(int run, double **X)
{
	int i;
	for(i = 0; i < NumVar; i++)
		SetBufValues(run, i, X[i]);
}

void S_System_Divide::SetBufValues(int run, int var, double* X)
{
	memcpy(bufValues[run].X[var], X, sizeof(double)*NumTrial );
}

void S_System_Divide::SetBufValues(int run, int var, int trial, double x)
{
	bufValues[run].X[var][trial] = x;
}
void S_System_Divide::NormalizeExpValues()
{
    int r, i, t;
    double min_val, max_val;
    double val, range;
    
    for(r = 0; r < NumRun+NumPerturb; r++)
    {
        for(i = 0; i < NumVar; i++)
        {
            //find min. and max. value
            min_val = expValues[r].X[i][0];
            max_val = min_val;
            
            for(t = 1; t < NumTrial; t++)
            {
                val	= expValues[r].X[i][t];
                if(val < min_val) min_val = val;
                if(val > max_val) max_val = val;
            }
            
            range = max_val - min_val;
            
            //normalization
            for(t = 0; t < NumTrial; t++)
            {
                val = expValues[r].X[i][t] - min_val;
                expValues[r].X[i][t] = (val <= 0 ? REPLACE_ZERO : val / range);
            }
        }//end of for each variable
    }//end of for each run
}
double S_System_Divide::ExpValue(int run, int i, int t) 
{
	return expValues[run].X[i][t];
}
double S_System_Divide::CalValue(int run, int i, int t) 
{
	return calValues[run].X[i][t];
}
double S_System_Divide::BufValue(int run, int i, int t)
{
	return bufValues[run].X[i][t];
}
double S_System_Divide::ExpDeltaValue(int run, int i, int t) 
{
	return expValues[run].dX[i][t];
}
double S_System_Divide::CalDeltaValue(int run, int i, int t) 
{
	return calValues[run].dX[i][t];
}
double S_System_Divide::CalFitnessValue(int run, int i, int t) 
{
	return calValues[run].fitness[i][t];
}
void S_System_Divide::SetGeneIndex(int index)
{
    m_GeneIndex=index;
}
