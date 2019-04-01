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
    best_cc = 0;
    function_cost = 0;
    g_NumConnections = 0;
    h_NumConnections = 0;
    time_interval = 0;
    min_fitness = XMAX;
    //printf("finish ini obj\n");
}
Objfunction::~Objfunction()
{
    int i;
    if(EVAL){delete [] EVAL;}
    if(GeneParamValues){delete [] GeneParamValues;}
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
    memset(EVAL,0,NO_POP);
    cal_profile = new double* [NumRuns];
    for(i=0;i<NumRuns;i++)
        cal_profile[i] = new double[NumTrials];
    ModelParameters = new double* [runtime];
    for(i=0;i<runtime;i++)
    {
        ModelParameters[i] = new double[NumGeneParam+1];
        memset(ModelParameters[i],0,NumGeneParam+1);
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
        if(pp>=0)
            GeneParamValues[1+pp]=genes[i+1];//Gij
    }
    for(i=2*g_NumConnections;i<(2*g_NumConnections+2*h_NumConnections-1);i+=2)//Hij
    {
        pp=(int)genes[i];
        //printf("pp=%d\n",pp);
        //PAUSE
        if(pp>=0)
            GeneParamValues[2+NumVars+pp]=genes[i+1];//Gij
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
    int j, r, i, x, ptr;
    double now_cc;
    double fitness, val,val2,val21,val22;
    double exp_avg,cal_avg,exp_total,cal_total,cal_diff,exp_diff,cal_diff_pow,exp_diff_pow,cal_exp_multi,c_c,c_c_fitness,add_t,exp_delta,interval_t;
    double sum, totalSum, ErrorValue;
    //double penalty;
    double *p,temp;
    double max_cc = 1;
    exp_diff_pow = 0;
    cal_diff_pow = 0;
    cal_exp_multi = 0;
    ptr = 0;
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
    for(r = 0; r < NumRuns; r++)
    {
        sum = 0;
        exp_total=0;
        cal_total=0;
        for(x = 0; x < NumTrials; x += time_interval)
        {
            exp_total += s_system->ExpValue(r, i, x);
            //printf("%f ",s_system->ExpValue(r, k, x));
        }
        for(x = 0; x < NumTrials; x += time_interval)
        {
            cal_total += s_system->CalValue(r, i, x);
            //printf("%f ",s_system->CalValue(r, k, x));
        }
        //printf("exp total:%f\n",exp_total);
        //printf("cal total:%f\n",cal_total);
        //PAUSE
        cal_exp_multi = 0;
        exp_diff_pow = 0;
        cal_diff_pow = 0;
        for(j = 0; j < NumTrials; j += time_interval)
        {
            val = s_system->ExpDeltaValue(r, i, j);//bechmark value
            //printf("exp_delta=%f\n",val);
            val2= s_system->CalDeltaValue(r, i, j);
            //printf("cal_delta=%f\n",val2);
            val21 = s_system->ExpValue(r, i, j);		//bechmark value
            add_t = s_system->ExpValue(r, i, j+1);
            exp_delta = (add_t-val21)/interval_t;
            //printf("cal_exp_delta=%f\n",exp_delta);
            //PAUSE
            val22= s_system->CalValue(r, i, j);
            //--------cc--------------------------//
            exp_avg = exp_total/NumTrials;
            cal_avg = cal_total/NumTrials;
            exp_diff = val21 - exp_avg;
            cal_diff = val22 - cal_avg;
            cal_exp_multi += exp_diff*cal_diff;
            exp_diff_pow += pow(exp_diff,2);
            cal_diff_pow += pow(cal_diff,2);
            //--------cc---------------------------//
            ErrorValue=fabs(val21 - val22);
            switch(fitness_type)
            {
                case 0:
                    sum += pow(1 - val22/val21,2);
                    break;
                case 1:
                    sum += pow(ErrorValue,2);
                    break;
                case 2:
                    temp=fabs(exp_delta - val2);
                    sum += (pow(temp,2) + pow(1 - val22/val21,2));
                    break;
                case 3:
                    temp=fabs(exp_delta - val2);
                    sum += (pow(temp,2) + pow(ErrorValue,2));
                    break;
                case 4:
                    temp=fabs(exp_delta - val2);
                    sum += (pow(temp,2));
                    break;
            }
        }//end of Numtrial
        totalSum += sum;
    }//end of NumRuns
    
    c_c = cal_exp_multi/(sqrt(cal_diff_pow)*sqrt(exp_diff_pow));
    //printf("cc:%f\n",c_c);
    if(c_c > 1){c_c = 1;}
    c_c_fitness = max_cc - c_c;
    totalSum /= (double)NumRuns;
    if(cc_mode)
    {
        fitness = totalSum + c_c_fitness;
    }
    else{
        fitness = totalSum;
    }
    now_cc = c_c;
    function_cost ++;
    if(fitness<min_fitness)
    {
        for(r = 0; r < NumRuns; r++)
        {
            for(j = 0; j < NumTrials; j += time_interval)
            {
                
                cal_profile[r][j]=s_system->CalValue(r, i, j);
                /*
                 switch(FUNC)
                 {
                 case 0:
                 tmp_sum += pow(1-cal_profile[j]/exp_profile[j],2);
                 break;
                 case 1:
                 tmp_sum += pow(cal_profile[j] - exp_profile[j],2);
                 break;
                 }
                 */
            }
        }//end of NumRuns
        min_fitness = fitness;
        //printf("min fitness=%e\n",min_fitness);
    }

    if(now_cc>=best_cc){best_cc = now_cc;}
    
    return fitness;
}

