#include "Objfunction.h"
extern int cc_mode;
Objfunction::Objfunction(S_System_Divide *ssystem,Knowledge *know)
{
    s_system = ssystem;
    knowledge = know;
    EVAL = NULL;
    cal_profile = NULL;
    GeneParamValues = NULL;
    ModelParameters = NULL;
    runtime =0;
    NumVars = 0;
    NumRuns = 0;
    NumTrials = 0;
    NumTFs = 0;
    PARAM_NUM = 0;
    NumGeneParam = 0;
    NO_POP = 0;
    last_cc = 0;
    RunCC = NULL;
    fitness_type = 2;
    function_cost = 0;
    g_NumConnections = 0;
    h_NumConnections = 0;
    time_interval = 0;
    //printf("finish ini obj\n");
}
Objfunction::~Objfunction()
{
    int i;
    if(EVAL){delete [] EVAL;}
    if(GeneParamValues){delete [] GeneParamValues;}
    if(RunCC){delete [] RunCC;}
    if(cal_profile)
    {
        for(i=0;i<NumRuns;i++)
            delete [] cal_profile[i];
        delete [] cal_profile;
    }
    if(ModelParameters)
    {
        for(i=0;i<runtime;i++)
            delete [] ModelParameters[i];
        delete [] ModelParameters;
    }

}
void Objfunction::init_obj(int time_step,int type,int run)
{
    int i;
    runtime = run;
    GeneParamValues = new double[NumGeneParam];
    EVAL = new double[NO_POP];
    memset(EVAL,0,NO_POP*sizeof(double));
    memset(GeneParamValues,0,NumGeneParam*sizeof(double));
    cal_profile = new double* [NumRuns];
    for(i=0;i<NumRuns;i++)
        cal_profile[i] = new double[NumTrials]();
    RunCC = new double[run]();
    ModelParameters = new double* [runtime];
    for(i=0;i<runtime;i++)
    {
        ModelParameters[i] = new double[NumGeneParam+1];
        memset(ModelParameters[i],0,(NumGeneParam+1)*sizeof(double));
    }
    g_NumConnections = knowledge->get_gNumconnections();
    h_NumConnections = knowledge->get_hNumconnections();
    time_interval = time_step;
    fitness_type = type;
}
void Objfunction::printCurrentModelParameters(double *genes,int type)
{
    int j,ii,ss;
    if(type==0)
        printf(".......................Decode \n");
    ii=0;
    ss=0;
    
    if(type==0)
    {
        for(j=0;j<PARAM_NUM;j++)
        {
            printf("%f, ",genes[ii+j]);
        }
        printf("\n");
    }
    for(j=0;j<NumGeneParam;j++)
    {
        printf("%f, ",GeneParamValues[ss+j]);
    }
    printf("\n");
    
    
}
void Objfunction::Decode(double *genes)
{
    int i;
    int pp;
    for(i=0;i<NumGeneParam;i++)
        GeneParamValues[i]=0;//orig
    
    GeneParamValues[0]=genes[PARAM_NUM-3];//alpha
    GeneParamValues[1+NumVars]=genes[PARAM_NUM-2];//beta
    GeneParamValues[NumGeneParam-1]=genes[PARAM_NUM-1];//interval_time
    
    for(i=0;i<2*g_NumConnections-1;i+=2)//Gij
    {
        pp=(int)genes[i];
        //printf("pp=%d\n",pp);
        //PAUSE
        if(pp>=0 && pp<NumVars)
            GeneParamValues[1+pp]=genes[i+1];//Gij
    }
    for(i=2*g_NumConnections;i<(2*g_NumConnections+2*h_NumConnections-1);i+=2)//Hij
    {
        pp=(int)genes[i];
        //printf("pp=%d\n",pp);
        //PAUSE
        if(pp>=0 && pp<NumVars)
            GeneParamValues[2+NumVars+pp]=genes[i+1];//Hij
    }
    
    //printCurrentModelParameters(genes,0);
    //PAUSE
}
double Objfunction::fitness(double *genes)
{
    double value;
    Decode(genes);
    //printCurrentModelParameters(genes,0);
    value=fitnessValue(GeneParamValues);
    return value;
}
double Objfunction::fitnessValue(double *genes)//(GeneParamValues)
{
    int j, r, i, x, n;
    double fitness, val,val2,val21,val22;
    double exp_avg,cal_avg,exp_total,cal_total,cal_diff,exp_diff,cal_diff_pow,exp_diff_pow,cal_exp_multi,c_c,c_c_fitness,add_t,exp_delta,interval_t;
    double sum, totalSum, ErrorValue, den;
    //double penalty;
    double *p,temp;
    double max_cc = 1;
    bool has_next;
    p=genes;
    //alpha
    i=s_system->m_GeneIndex;
    val=p[0];
    s_system->SetCalAlpha(i, val);
    //Gij
    for(j = 0; j < NumVars; j++)
    {
        val=p[j+1];//if NumVars=5, p[1],p[2],p[3],...p[5]
        s_system->SetCalG(i, j, val);//
    }
    //Beta
    val=p[1+NumVars];
    s_system->SetCalBeta(i,val );
    //Hij
    for(j = 0; j < NumVars; j++)
    {
        val=p[2+NumVars+j];//*p[3*(NumVars+1)+1+j];//p[j+1]+;
        s_system->SetCalH(i, j, val);
        
    }
    //delta[j]
    val=p[NumGeneParam-1];
    s_system->SetDeltaT(val);
    
    interval_t=s_system->GetDeltaT();
    fitness = 0;
    
    s_system->CalAllCalValues(); //calculate the estimate profile
    
    //Expression Error
    totalSum = 0;
    //cc: Pearson correlation of the observed and simulated profile over ALL replicates. Every
    //replicate is centred with its own mean (over the sampled time points), the sums are pooled.
    cal_exp_multi = 0;
    exp_diff_pow = 0;
    cal_diff_pow = 0;
    for(r = 0; r < NumRuns; r++)
    {
        sum = 0;
        exp_total=0;
        cal_total=0;
        n = 0;
        for(x = 0; x < NumTrials; x += time_interval)
        {
            exp_total += s_system->ExpValue(r, i, x);
            cal_total += s_system->CalValue(r, i, x);
            n++;
        }
        //the mean is taken over the sampled points (the sums only contain every time_interval-th point)
        exp_avg = exp_total/n;
        cal_avg = cal_total/n;
        for(j = 0; j < NumTrials; j += time_interval)
        {
            val = s_system->ExpDeltaValue(r, i, j);//bechmark value
            val2= s_system->CalDeltaValue(r, i, j);
            val21 = s_system->ExpValue(r, i, j);		//bechmark value
            //The observed derivative needs the next time point. For observed (not simulated) data
            //there is none after the last point (X[NumTrials] is only allocated, never measured),
            //so the derivative term is skipped there.
            has_next = (j+1 < NumTrials) || s_system->IsSimulated();
            if(has_next)
            {
                add_t = s_system->ExpValue(r, i, j+1);
                exp_delta = (add_t-val21)/interval_t;
            }
            else
                exp_delta = 0;
            val22= s_system->CalValue(r, i, j);
            //--------cc--------------------------//
            exp_diff = val21 - exp_avg;
            cal_diff = val22 - cal_avg;
            cal_exp_multi += exp_diff*cal_diff;
            exp_diff_pow += exp_diff*exp_diff;
            cal_diff_pow += cal_diff*cal_diff;
            //--------cc---------------------------//
            ErrorValue=fabs(val21 - val22);
            temp = has_next ? fabs(exp_delta - val2) : 0;
            switch(fitness_type)
            {
                case 0:
                    sum += pow(1 - val22/val21,2);
                    break;
                case 1:
                    sum += pow(ErrorValue,2);
                    break;
                case 2:
                    sum += (pow(temp,2) + pow(1 - val22/val21,2));
                    break;
                case 3:
                    sum += (pow(temp,2) + pow(ErrorValue,2));
                    break;
                case 4:
                    sum += (pow(temp,2));
                    break;
            }
        }//end of Numtrial
        totalSum += sum;
    }//end of NumRuns
    
    //a constant observed or simulated profile has no correlation: cc = 0 (it was 0/0 = NaN)
    den = sqrt(cal_diff_pow)*sqrt(exp_diff_pow);
    if(den > 0 && den < HUGE_VAL)
        c_c = cal_exp_multi/den;
    else
        c_c = 0;
    if(!(c_c == c_c)){c_c = 0;}
    if(c_c > 1){c_c = 1;}
    if(c_c < -1){c_c = -1;}
    c_c_fitness = max_cc - c_c;
    totalSum /= (double)NumRuns;
    if(cc_mode)
    {
        fitness = totalSum + c_c_fitness;
    }
    else{
        fitness = totalSum;
    }
    //NaN is the worst fitness (it would otherwise never compare as worse and could stay "best")
    if(!(fitness == fitness)){fitness = NAN_PENALTY;}
    last_cc = c_c;
    function_cost ++;
    return fitness;
}
void Objfunction::RecordRun(int run,double *genes)
{
    int r, j;
    unsigned int cost = function_cost;
    fitness(genes);//the simulated profile and cc of this chromosome
    function_cost = cost;
    RunCC[run] = last_cc;
    for(r = 0; r < NumRuns; r++)
    {
        for(j = 0; j < NumTrials; j += time_interval)
            cal_profile[r][j]=s_system->CalValue(r, s_system->m_GeneIndex, j);
    }
}
