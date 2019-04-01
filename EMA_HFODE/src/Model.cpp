#include "Model.h"

Model::Model()//initialization
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
    //m_GeneIndex = 0;
    //CheckIndex = 0;
	//experiment
    expb = NULL;
    expTransMax = NULL;
    //exp_beta = NULL;
    exp_n = NULL;
    exp_k = NULL;
    expdegrade = NULL;
	isSetInitExpVal = NULL;
	expValues = NULL;
    isInitExpModel = false;
	//calculated
    calb = NULL;
    calTransMax = NULL;
    //cal_beta = NULL;
    cal_n = NULL;
    cal_k = NULL;
    caldegrade = NULL;
    isSetInitCalVal = NULL;
    calValues = NULL;
    isInitCalModel = false;
    //range
    //deltat_up = 0;
    //deltat_down = 0;
    b_up = 0;
    b_down = 0;
    TransMax_up = 0;
    TransMax_down = 0;
    //beta_up = 0;
    //beta_down = 0;
    n_up = 0;
    n_down = 0;
    k_up = 0;
    k_down = 0;
    degrade_up = 0;
    degrade_down = 0;
    count = 0;
    //printf("finish init s\n");
}
Model::~Model()
{
	int i, j;
    if(NumVar != 0)
    {
        //experiment
        if(expb) delete [] expb;
        if(expdegrade)  delete [] expdegrade;
        if(expTransMax) delete [] expTransMax;
        /*
        if(exp_beta)
        {
            for(i = 0; i < NumVar; i++)
                delete [] exp_beta[i];
            delete [] exp_beta;
        }
         */
        if(exp_n)
        {
            for(i = 0; i < NumVar; i++)
                delete [] exp_n[i];
            delete [] exp_n;
        }
        if(exp_k)
        {
            for(i = 0; i < NumVar; i++)
                delete [] exp_k[i];
            delete [] exp_k;
        }
        //calculated
        if(calb) delete [] calb;
        if(caldegrade)  delete [] caldegrade;
        if(calTransMax) delete [] calTransMax;
        /*
        if(cal_beta)
        {
            for(i = 0; i < NumVar; i++)
                delete [] cal_beta[i];
            delete [] cal_beta;
        }
         */
        if(cal_n)
        {
            for(i = 0; i < NumVar; i++)
                delete [] cal_n[i];
            delete [] cal_n;
        }
        if(cal_k)
        {
            for(i = 0; i < NumVar; i++)
                delete [] cal_k[i];
            delete [] cal_k;
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
                //delete [] expValues[i].fitness[j];
            }
            delete [] expValues[i].X;
            delete [] expValues[i].dX;
            //delete [] expValues[i].fitness;
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
                //delete [] calValues[i].fitness[j];
			}
			delete [] calValues[i].X;
			delete [] calValues[i].dX;
            //delete [] calValues[i].fitness[j];
		}
		delete [] calValues;
	}

	//Experiment
	if(isSetInitExpVal != NULL)
		delete [] isSetInitExpVal;
    
	//calculated
	if(isSetInitCalVal != NULL)
		delete [] isSetInitCalVal;
    
}

void Model::InitializeExp()
{
    //printf("InitializeExp\n");
    int i, j;
    expb = new double[NumVar];
    memset(expb, 0, NumVar*sizeof(double));
    expdegrade = new double[NumVar];
    memset(expdegrade, 0, NumVar*sizeof(double));
    expTransMax = new double[NumVar];
    memset(expTransMax, 0, NumVar*sizeof(double));
    //exp_beta = new double* [NumVar];
    //exp_n = new int* [NumVar];
    exp_n = new double* [NumVar];
    exp_k = new double* [NumVar];
    for(i = 0; i < NumVar; i++)
    {
        //exp_beta[i] = new double[NumVar];
        //memset(exp_beta[i], 0, NumVar*sizeof(double));
        //exp_n[i] = new int[NumVar];
        exp_n[i] = new double[NumVar];
        memset(exp_n[i], 0, NumVar*sizeof(double));
        exp_k[i] = new double[NumVar];
        memset(exp_k[i], 0, NumVar*sizeof(double));
    }
    expValues = new Values[NumRun+NumPerturb];
    //printf("half\n");
    for(i = 0; i < NumRun+NumPerturb; i++)
    {
        expValues[i].X = new double* [NumVar];
        expValues[i].dX = new double* [NumVar];
        //expValues[i].fitness = new double* [NumVar];
        
        for(j = 0; j < NumVar; j++)
        {
            expValues[i].X[j] = new double[NumTrial+1];
            expValues[i].dX[j] = new double[NumTrial+1];
            //expValues[i].fitness[j] = new double[NumTrial+1];
        }
    }
    isSetInitExpVal = new bool[NumRun+NumPerturb];
    for(i = 0; i < NumRun+NumPerturb; i++)
        isSetInitExpVal[i] = false;
    
    isInitExpModel = true;
}

void Model::InitializeCal()
{
    //printf("InitializeCal\n");
    int i, j;
    calb = new double[NumVar];
    memset(calb, 0, NumVar*sizeof(double));
    caldegrade = new double[NumVar];
    memset(caldegrade, 0, NumVar*sizeof(double));
    calTransMax = new double[NumVar];
    memset(calTransMax, 0, NumVar*sizeof(double));
    //cal_beta = new double* [NumVar];
    //cal_n = new int* [NumVar];
    cal_n = new double* [NumVar];
    cal_k = new double* [NumVar];
    for(i = 0; i < NumVar; i++)
    {
        //cal_beta[i] = new double[NumVar];
        //memset(cal_beta[i], 0, NumVar*sizeof(double));
        //cal_n[i] = new int[NumVar];
        cal_n[i] = new double[NumVar];
        memset(cal_n[i], 0, NumVar*sizeof(double));
        cal_k[i] = new double[NumVar];
        memset(cal_k[i], 0, NumVar*sizeof(double));
    }
    calValues = new Values[NumRun+NumPerturb];
    for(i = 0; i < NumRun+NumPerturb; i++)
    {
        calValues[i].X = new double* [NumVar];
        calValues[i].dX = new double* [NumVar];
        //calValues[i].fitness = new double* [NumVar];
        
        for(j = 0; j < NumVar; j++)
        {
            calValues[i].X[j] = new double[NumTrial+1];
            calValues[i].dX[j] = new double[NumTrial+1];
            //calValues[i].fitness[j] = new double[NumTrial+1];
        }
    }
    isSetInitCalVal = new bool[NumRun+NumPerturb];
    for(i = 0; i < NumRun+NumPerturb; i++)
        isSetInitCalVal[i] = false;
    
    isInitCalModel = true;
}
double Model::AddGaussianNoise(double stdev, double mean)
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

int Model::ReadInitFromFile(char* filename)
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
		else if(strcmp(token, "delta_t") == 0)
		{
			token = strtok(NULL, "\t\r\n #");
			delta_t = atof(token);
            //printf("delta_t=%f\n",delta_t);
            //PAUSE
		}
        /*
        else if(strcmp(token, "delta_s") == 0)
		{
			token = strtok(NULL, "\t\r\n #");
			delta_s = atof(token);
            //printf("delta_s=%f\n",delta_s);
		}
         */
        else if(strcmp(token, "nochange") == 0)
        {
            token = strtok(NULL, "\t\r\n #");
            Nochange = atoi(token);
            //printf("Nochange=%d\n",Nochange);
        }
        else if(strcmp(token, "b_up") == 0)
        {
            token = strtok(NULL, "\t\r\n #");
            b_up = atoi(token);
            //printf("b_up=%f\n",b_up);
        }
        else if(strcmp(token, "b_down") == 0)
        {
            token = strtok(NULL, "\t\r\n #");
            b_down = atoi(token);
            //printf("b_down=%f\n",b_down);
        }
        else if(strcmp(token, "transMax_up") == 0)
        {
            token = strtok(NULL, "\t\r\n #");
            TransMax_up = atoi(token);
        }
        else if(strcmp(token, "transMax_down") == 0)
        {
            token = strtok(NULL, "\t\r\n #");
            TransMax_down = atoi(token);
        }
        /*
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
         */
        else if(strcmp(token, "n_up") == 0)
        {
            token = strtok(NULL, "\t\r\n #");
            n_up = atoi(token);
            //printf("n_up=%d\n",n_up);
        }
        else if(strcmp(token, "n_down") == 0)
        {
            token = strtok(NULL, "\t\r\n #");
            n_down = atoi(token);
            //printf("n_downn=%d\n",n_down);
        }
        else if(strcmp(token, "k_up") == 0)
        {
            token = strtok(NULL, "\t\r\n #");
            k_up = atoi(token);
            //printf("k_up=%f\n",k_up);
        }
        else if(strcmp(token, "k_down") == 0)
        {
            token = strtok(NULL, "\t\r\n #");
            k_down = atoi(token);
            //printf("k_down=%f\n",k_down);
        }
        else if(strcmp(token, "degrade_up") == 0)
        {
            token = strtok(NULL, "\t\r\n #");
            degrade_up = atoi(token);
            //printf("degrade_up=%f\n",degrade_up);
        }
        else if(strcmp(token, "degrade_down") == 0)
        {
            token = strtok(NULL, "\t\r\n #");
            degrade_down = atoi(token);
            //printf("degrade_down=%f\n",degrade_down);
            //PAUSE
        }
        /*
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
         */
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
                //b
                token = strtok(line, "\n\r\t #");
                expb[i] = atof(token);
                //printf("expb[%d]=%f\n",i,expb[i]);
                token = strtok(NULL, "\n\r\t #");
                expTransMax[i] = atof(token);
                //printf("TransMax[%d] = %f\n",i, expTransMax[i]);
                
                //Nij
                for(j = 0; j < NumVar; j++)
                {
                    token = strtok(NULL, "\n\r\t #");
                    exp_n[i][j] = atoi(token);
                    //printf("exp_n[%d][%d]=%f\n",i,j,exp_n[i][j]);
                }
                
                /*
                //betaij
                for(j = 0; j < NumVar; j++)
                {
                    token = strtok(NULL, "\n\r\t #");
                    exp_beta[i][j] = atoi(token);
                    //printf("exp_g[%d][%d]=%f\n",i,j,exp_g[i][j]);
                }
                */
                //Kij
                for(j = 0; j < NumVar; j++)
                {
                    token = strtok(NULL, "\n\r\t #");
                    exp_k[i][j] = atof(token);
                    //printf("exp_k[%d][%d]=%f\n",i,j,exp_k[i][j]);
                }
                //degrade
                token = strtok(NULL, "\n\r\t #");
                expdegrade[i] = atof(token);
                //printf("expdegrade[%d]=%f\n",i,expdegrade[i]);
            }
			isSimulated = true;
            //printf("finish to readd the init parameter\n");
		}
        
        
		else if(strcmp(token, "list_data") == 0)
		{
            //printf("list_data\n");
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
void Model::PrintAllExpXs(FILE* fout)
{
	int r;
	int cnt = NumRun + NumPerturb;
    //printf("cnt=%d\n",cnt);
	for(r = 0; r < cnt; r++)
	{
		PrintExpXs(fout, r);
        
		fprintf(fout, "\n\n");
	}
}

void Model::PrintExpXs(FILE* fout, int run)
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
        //printf("\n");
		fprintf(fout, "\n");
	}
}

void Model::PrintExpXs(FILE* fout, int run, int var)
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

void Model::PrintExpX(FILE* fout, int run, int var, int trial)
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

void Model::SetInitExpAllRun(double** X) //X[r][i][0]: run r, var i, trial 0
{
	int i;
	for(i = 0; i < NumRun; i++)
		SetInitExpRun(i, X[i]);
}

void Model::SetInitExpRun(int nRun, double* X)
{
	int i;
	for(i = 0; i < NumVar; i++)
		expValues[nRun].X[i][0] = X[i];
	
	isSetInitExpVal[nRun] = true;
}

void Model::CalAllExpValues()//simulated exp
{
	//if(!isSimulated)
	//	return -1;
    
	int i, v, t, j;
	double tmp;
	//int rtn = 1;
	for(i = 0; i < NumRun; i++)
	{
		//rtn *= CalRunExpValues(i);
        CalRunExpValues(i);
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
    
	//return rtn;
}

void Model::CalRunExpValues(int nRun)//produce simulated profile
{
    double tmp, temp, activate;
    int i, j, t;
    //double expval, val, X;
    
    if(isSetInitExpVal[nRun] != true)
    {
        fprintf(stderr, "Experiment Run %d is not initialized!\n", nRun);
        exit(0);
    }
    
    
    for(t = 0; t < NumTrial; t++)
    {
        for(i = 0; i < NumVar; i++)
        {
            temp=1;
            //activate=0;
            for(j = 0; j < NumVar; j++)
            {
                //printf("exp_n[%d][%d]=%f\texp_k[%d][%d]=%f\texpValue[%d].[%d][%d]=%f",i,j,exp_n[i][j],i,j,exp_k[i][j],nRun,j,t,expValues[nRun].X[j][t]);PAUSE;
                if(exp_n[i][j] == 0){
                    temp *= 1;
                    //printf("temp=%f\t",temp);PAUSE;
                }else{
                    temp *= 1/(1+pow(exp_k[i][j]/expValues[nRun].X[j][t],exp_n[i][j]));
                    //printf("temp=%f\t",temp);PAUSE;
                }
                    
            }
            activate = expTransMax[i]*temp + expb[i];
            //printf("active=%f\n",activate);
            if(activate<0)
            {
                activate = 0;
            }
            expValues[nRun].dX[i][t] = activate - expdegrade[i]*expValues[nRun].X[i][t];
            tmp = expValues[nRun].dX[i][t] * delta_t + expValues[nRun].X[i][t];
            //printf("tmp=%f\n",tmp);
            if(tmp >= XMAX)
                expValues[nRun].X[i][t+1] = XMAX;
            else if(tmp <= XMIN)
                expValues[nRun].X[i][t+1] = XMIN;
            else
                expValues[nRun].X[i][t+1] = tmp;
            //printf("exp=%f\n",expValues[nRun].X[i][t+1]);
        }//end of NumVar
    }//end of number of trial
}

void Model::SetExpValues(int r, int i, int t, double val)
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

void Model::SetExpb(int i, double val)
{
    if(i < NumVar && expb != NULL)
        expb[i] = val;
    else
        fprintf(stderr, "Error: SetExpb out of index range!!!\n");
}

void Model::SetExpdegrade(int i, double val)
{
    if(i < NumVar && expdegrade != NULL)
        expdegrade[i] = val;
    else
        fprintf(stderr, "Error: SetExpdegrade out of index range!!!\n");
}
void Model::SetExpTransMax(int i , double val)
{
    if(i < NumVar && expTransMax != NULL)
        expTransMax[i] = val;
    else
        fprintf(stderr, "Error: SetExpTransMax out of index range!!!\n");
}
/*
void Model::SetExpbeta(int i, int j, double val)
{
    if(i < NumVar && j < NumVar && exp_beta != NULL)
        exp_beta[i][j] = val;
    else
        fprintf(stderr, "Error: SetExpbeta out of index range!!!\n");
}
*/
void Model::SetExpn(int i, int j, int val)
{
    if(i < NumVar && j < NumVar && exp_n != NULL)
        exp_n[i][j] = val;
    else
        fprintf(stderr, "Error: SetExpn out of index range!!!\n");
}

void Model::SetExpk(int i, int j, double val)
{
    if(i < NumVar && j < NumVar && exp_k != NULL)
        exp_k[i][j] = val;
    else
        fprintf(stderr, "Error: SetExp out of index range!!!\n");
}

void Model::SetInitCalAllRun(double** X) //X[r][i][0]: run r, var i, trial 0
{
	int i;
	for(i = 0; i < NumRun; i++)
		SetInitCalRun(i, X[i]);
}

void Model::SetInitCalRun(int nRun, double* X)
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

void Model::PrintAllCalXs(FILE* fout)
{
	int r;
	int cnt = ExecutionRuns;
	
	for(r = 0; r < cnt; r++)
	{
		PrintCalXs(fout, r);
        //fprintf(fout, "\n\n");
	}
}

void Model::PrintCalXs(FILE* fout, int run)
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

void Model::PrintCalXs(FILE* fout, int run, int var)
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

void Model::PrintCalX(FILE* fout, int run, int var, int trial)
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
*/
void Model::CalAllCalValues(int index)//calculate the fitting profile
{
	int i, v;
	//int rtn = 1;
    int totalnumRun = NumRun + NumPerturb;
    //printf("totalnumRun=%d\n",totalnumRun);
    //PAUSE
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
        //rtn *=CalRunVarCalValues_Divide(i, index);//, j);	//run times, number of xi
        CalRunVarCalValues_Divide(i, index);
	}
    //return rtn;
}
void Model::CalRunVarCalValues_Divide(int nRun, int GeneIndex)// calculate the fitting profile
{
    double temp, tmp, activate, tempx;
    int i, j, t;
    //tmp = 0;
    //tempx = 0;
    //printf("m_GeneIndex=%d\tCheckIndex=%dGeneIndex=%d\n",m_GeneIndex,CheckIndex,GeneIndex);
    //if(CheckIndex != m_GeneIndex)
    //    m_GeneIndex = CheckIndex;
    //int flag;
    //int const Var = m_GeneIndex;
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
        exit(0);
    }
    
    
    for(t = 0; t < NumTrial; t++)//Time points
    {
        //Var=m_GeneIndex;
        
        temp=1;
        //temp = 1;
        //activate = 0;
        for(j = 0; j < NumVar; j++)	//
        {
            //printf("cal_n[%d][%d]:%f\t",GeneIndex,j,cal_n[GeneIndex][j]);
            //printf("cal_k[%d][%d]:%f\n",GeneIndex,j,cal_k[GeneIndex][j]);
            //PAUSE
            //Divide and conquer technique
            if(j == GeneIndex)
            {
                tempx= calValues[nRun].X[j][t];
            }
            else//if(i != j)
            {
                tempx = expValues[nRun].X[j][t];
                //printf("tempx=%f\tnRun=%d\tj=%d\tt=%d\n",tempx,nRun,j,t);
            }
            if(cal_n[GeneIndex][j] == 0)
            {
                temp *= 1;
            }else{
                temp *= 1/(1+pow(cal_k[GeneIndex][j]/tempx,cal_n[GeneIndex][j]));
            }
            //temp += cal_beta[GeneIndex][j]/(1+pow(cal_k[GeneIndex][j]/tempx,cal_n[GeneIndex][j]));
            //temp *= 1+cal_beta[Var][j]/(1+pow(cal_k[Var][j]/tempx,cal_n[Var][j]));
            
        }
        //printf("temp=%f\n",temp);
        activate = calTransMax[GeneIndex]*temp + calb[GeneIndex];
        //printf("activate:%f\n",activate);
        //PAUSE
        if(activate<0){
            activate = 0;
        }
        calValues[nRun].dX[GeneIndex][t] = activate - caldegrade[GeneIndex] * calValues[nRun].X[GeneIndex][t];
        tmp = calValues[nRun].dX[GeneIndex][t] * delta_t + calValues[nRun].X[GeneIndex][t];

        //printf("%d tmp=%f\n",t,tmp);
        if(tmp >= XMAX)
            calValues[nRun].X[GeneIndex][t+1] = XMAX;
        else if(tmp <= XMIN)
            calValues[nRun].X[GeneIndex][t+1] = XMIN;
        else
            calValues[nRun].X[GeneIndex][t+1] = tmp;
    }//end of number of trial
    //return 1;
    count++;
}
void Model::SetCalb(int i, double val)
{
    if(i < NumVar && calb != NULL)
        calb[i] = val;
    else
    {
       fprintf(stderr, "i:%d, Error: SetCalb out of index range!!!\n",i);
        exit(0);
    }
    
}

void Model::SetCaldegrade(int i, double val)
{
    if(i < NumVar && caldegrade != NULL)
        caldegrade[i] = val;
    else
    {
        fprintf(stderr, "i:%d, Error: SetCaldegrade out of index range!!!\n",i);
        exit(0);
    }
    
}
void Model::SetCalTransMax(int i, double val)
{
    if(i < NumVar && calTransMax != NULL)
        calTransMax[i] = val;
    else
    {
        fprintf(stderr, "i:%d, Error: SetCalTransMax out of index range!!!\n",i);
        exit(0);
    }
}
/*
void Model::SetCalbeta(int i, int j, double val)
{
    if(i < NumVar && j < NumVar && cal_beta != NULL)
        cal_beta[i][j] = val;
    else
    {
        fprintf(stderr, "i:%d,j:%d, Error: SetCalbeta out of index range!!!\n",i,j);
        exit(0);
    }
    
}
*/
void Model::SetCaln(int i, int j, double val)
{
    if(i < NumVar && j < NumVar && cal_n != NULL)
        cal_n[i][j] = val;
    else
    {
        fprintf(stderr, "i:%d,j:%d, Error: SetCaln out of index range!!!\n",i,j);
        exit(0);
    }
    
}
void Model::SetCalk(int i, int j, double val)
{
    if(i < NumVar && j < NumVar && cal_k != NULL)
        cal_k[i][j] = val;
    else
    {
        fprintf(stderr, "i:%d,j:%d,Error: SetCalk out of index range!!!\n",i,j);
        exit(0);
    }
        
    
}


void Model::NormalizeExpValues()
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
double Model::ExpValue(int run, int i, int t)
{
	return expValues[run].X[i][t];
}
double Model::CalValue(int run, int i, int t)
{
	return calValues[run].X[i][t];
}
//double Model::ExpDeltaValue(int run, int i, int t)
//{
//	return expValues[run].dX[i][t];
//}
double Model::CalDeltaValue(int run, int i, int t)
{
	return calValues[run].dX[i][t];
}
//double S_System_Divide::CalFitnessValue(int run, int i, int t)
//{
//	return calValues[run].fitness[i][t];
//}
//void Model::SetGeneIndex(int index)
//{
//    m_GeneIndex = index;
//}
