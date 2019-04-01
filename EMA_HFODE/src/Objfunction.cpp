#include "Objfunction.h"
//extern int cc_mode;
//Objfunction::Objfunction(Model *hill_f)
Objfunction::Objfunction()
{
    //hill = hill_f;
    //hill = new Model();
    //knowledge = new Getknowledge();
    //EVAL = NULL;
    cal_profile = NULL;
    GeneParamValues = NULL;
    Chromosomes = NULL;
    hill = NULL;
    //ModelParameters = NULL;
    NumVars = 0;
    NumRuns = 0;
    NumTrials = 0;
    NumTFs = 0;
    PARAM_NUM = 0;
    NumGeneParam = 0;
    NO_POP = 0;
    best_cc = 0;
    function_cost = 0;
    NumConnections = 0;
    //h_NumConnections = 0;
    time_interval = 0;
    runtime = 0;
    min_fitness = XMAX;
    //printf("finish ini obj\n");
}
Objfunction::~Objfunction()
{
    int i;
    //if(EVAL){delete [] EVAL;}
    if(GeneParamValues){delete [] GeneParamValues;}
    if(cal_profile)
    {
        for(i=0;i<NumRuns;i++)
            delete [] cal_profile[i];
        delete [] cal_profile;
    }
    /*
    if(hill)
    {
        delete hill;
    }
    if(Chromosomes)
    {
        delete Chromosomes;
    }
    */
}
void Objfunction::init_obj(int time_step,int type,int run, int cc, int No_reg, int index, Model *hill_f, double **Chrom)
{
    hill = hill_f;
    int i;
    runtime = run;
    cc_mode = cc;
    Chromosomes = Chrom;
    GeneParamValues = new double[NumGeneParam];
    //EVAL = new double[NO_POP];
    //memset(EVAL,0,NO_POP*sizeof(double));
    cal_profile = new double* [NumRuns];
    for(i=0;i<NumRuns;i++)
        cal_profile[i] = new double[NumTrials];
    NumConnections = No_reg;
    time_interval = time_step;
    fitness_type = type;
    GeneIndex = index;
    //ModelParameters = new double* [runtime];
    //for(i=0;i<runtime;i++)
    //{
    //  ModelParameters[i] = new double[NumGeneParam+1];
    //    memset(ModelParameters[i],0,NumGeneParam+1);
    //}

}
void Objfunction::printCurrentModelParameters(int PopIndex)
{
    int j;
    //ii,ss;
    printf("===============Decode start==========\n");
    printf("===========Encode Chromosome=========\n");
    //ii=0;
    //ss=0;
    for(j=0;j<PARAM_NUM;j++)
    {
        printf("%f, ",Chromosomes[PopIndex][j]);
    }
    printf("\n");
    //Decode(genes, PopIndex);
    printf("==========Decode Chromosome==========\n");
    for(j=0; j<NumGeneParam; j++)
    {
        printf("%f, ",GeneParamValues[j]);
    }
    printf("\n");
}
void Objfunction::Decode(int PopIndex)
{
    int i,j;
    int pp,index;
    int MM, mask;
    //printf("PopIndex=%d\n",PopIndex);
    //printf("bi=%f\n",chrom[PopIndex][PARAM_NUM-4]);
    //printf("degrade=%f\n",chrom[PopIndex][PARAM_NUM-3]);
    //printf("Mask=%f\n",chrom[PopIndex][PARAM_NUM-2]);
    for(i=0;i<NumGeneParam;i++){
        GeneParamValues[i] = 0;//initial
    }
    //GeneParamValues is origin chromosome, Chromosome is the hybrid encode, MM is mask
    GeneParamValues[0] = Chromosomes[PopIndex][PARAM_NUM-4];//bi
    GeneParamValues[1] = Chromosomes[PopIndex][PARAM_NUM-3];//TransMax
    GeneParamValues[NumGeneParam-1] = Chromosomes[PopIndex][PARAM_NUM-2];//degrade
    //GeneParamValues[NumGeneParam-1] = Chromosomes[PopIndex][PARAM_NUM-1];//interval_time
    MM = (int)Chromosomes[PopIndex][PARAM_NUM-1];//Mask
    j=0;
    index=0;
    for(i=0;i<2*NumConnections-1;i+=2,j++)
    {
        pp = (int)Chromosomes[PopIndex][i];
        mask = 1<<j;
        //printf("pp=%d\n",pp);
        //PAUSE
        if(pp >= 0)
        {
            if(MM & mask)
            {
                GeneParamValues[2+pp] = Chromosomes[PopIndex][i+1];//n
                GeneParamValues[2+NumVars+pp] = Chromosomes[PopIndex][2*NumConnections+index];//k
               // GeneParamValues[1+2*NumVars+pp] = Chromosomes[PopIndex][3*NumConnections+index];//k
            }
        }else{
            if(pp == neg_zero)
                pp = 0;
            else
                pp = abs(pp);
            GeneParamValues[2+pp] = 0;
            GeneParamValues[2+NumVars+pp] = 0;
            //GeneParamValues[1+2*NumVars+pp] = 0;
        }
        index++;
    }
    //printCurrentModelParameters(PopIndex);
    //PAUSE;
}
double Objfunction::OAfitenss(double **OATable, int OATableIndex)
{
    int i,j;
    int pp,index;
    int MM, mask;
    double value = 0;
    for(i=0;i<NumGeneParam;i++){
        GeneParamValues[i] = 0;//initial
    }
    
    GeneParamValues[0] = OATable[OATableIndex][PARAM_NUM-4];//bi
    GeneParamValues[1] = OATable[OATableIndex][PARAM_NUM-3];//TransMax
    //GeneParamValues[NumGeneParam-2] = OATable[OATableIndex][PARAM_NUM-3];//degrade
    GeneParamValues[NumGeneParam-1] = OATable[OATableIndex][PARAM_NUM-2];//degrade
    MM = (int)OATable[OATableIndex][PARAM_NUM-1];//Mask
    j=0;
    index=0;
    for(i=0;i<2*NumConnections-1;i+=2,j++)
    {
        pp=(int)OATable[OATableIndex][i];
        mask = 1<<j;
        //printf("pp=%d\n",pp);
        //PAUSE
        if(pp >= 0)
        {
            if(MM & mask)
            {
                GeneParamValues[2+pp] = OATable[OATableIndex][i+1];//n
                GeneParamValues[2+NumVars+pp] = OATable[OATableIndex][2*NumConnections+index];//k
                //GeneParamValues[1+2*NumVars+pp] = OATable[OATableIndex][3*NumConnections+index];//k
            }
        }else{
            if(pp == neg_zero)
                pp = 0;
            else
                pp = abs(pp);
            GeneParamValues[2+pp] = 0;
            GeneParamValues[2+NumVars+pp] = 0;
            //GeneParamValues[1+2*NumVars+pp] = 0;
        }
        index++;
    }
    /*
    for(j=0; j<NumGeneParam; j++)
    {
        printf("%f, ",GeneParamValues[j]);
    }
    printf("\n");
    PAUSE;
    */
    value = fitnessValue(GeneParamValues);
    //printf("fitness=%e\n",value);
    return value;
}
double Objfunction::fitness(int PopIndex)
{
    double value = 0;
    Decode(PopIndex);
    //printCurrentModelParameters(genes,0);
    value=fitnessValue(GeneParamValues);
    return value;
}
double Objfunction::fitnessValue(double *OneChrom)//(GeneParamValues)
{
    int j, r, x;
    //ptr;
    //int i;
    double now_cc;
    double fitness, val,val2,val21,val22;
    double exp_avg,cal_avg,exp_total,cal_total,cal_diff,exp_diff,cal_diff_pow,exp_diff_pow,cal_exp_multi,c_c,c_c_fitness,add_t,exp_delta,interval_t;
    double sum, totalSum, ErrorValue;
    //double penalty;
    //double *p = NULL;
    double temp;
    double max_cc = 1;
    exp_diff_pow = 0;
    cal_diff_pow = 0;
    cal_exp_multi = 0;
    //i = GeneIndex;
    //printf("Numvar=%d\n",NumVars);
    //printf("i=%d\n",i);
    //PAUSE
    val = OneChrom[0];
    //b
    //printf("b=%f\t",val);
    hill->SetCalb(GeneIndex, val);
    //TransMax
    val = OneChrom[1];
    //printf("TransMax=%f\t",val);
    hill->SetCalTransMax(GeneIndex, val);
    //Nij
    //printf("Nij:");
    for(j = 0; j < NumTFs; j++)
    {
        val = OneChrom[j+2];
        //printf("%f\t",val);
        hill->SetCaln(GeneIndex, j, val);
    }
    /*
    //Betaij
    for(j = 0; j < NumTFs; j++)
    {
        val=OneChrom[1+NumVars+j];
        hill->SetCalbeta(GeneIndex, j, val);
    }
     */
    //Kij
    //printf("kij:");
    for(j = 0; j < NumTFs; j++)
    {
        val = OneChrom[2+NumVars+j];
        //printf("%f\t",val);
        hill->SetCalk(GeneIndex, j, val);
    }
    //degrade[i]
    val = OneChrom[NumGeneParam-1];
    hill->SetCaldegrade(GeneIndex, val);
    //time interval[j]
    //val=OneChrom[NumGeneParam-1];
    //hill->SetDeltaT(val);
    interval_t=hill->GetDeltaT();
    //fitness = 0;
    
    hill->CalAllCalValues(GeneIndex); //calculate the estimate profile
    
    //Expression Error
    totalSum = 0;
    for(r = 0; r < NumRuns; r++)
    {
        sum = 0;
        
        
        //i=hill->m_GeneIndex;
        exp_total=0;
        cal_total=0;
        for(x = 0; x < NumTrials; x += time_interval)
            exp_total += hill->ExpValue(r, GeneIndex, x);
    
        for(x = 0; x < NumTrials; x += time_interval)
            cal_total += hill->CalValue(r, GeneIndex, x);
        
        //printf("exp total:%f\n",exp_total);
        //printf("cal total:%f\n",cal_total);
        //PAUSE
        cal_exp_multi = 0;
        exp_diff_pow = 0;
        cal_diff_pow = 0;
        for(j = 0; j < NumTrials; j += time_interval)
        {
            //val = s_system->ExpDeltaValue(r, k, j);//bechmark value
            //printf("exp_delta=%f\n",val);
            val2= hill->CalDeltaValue(r, GeneIndex, j);
            //printf("cal_delta=%f\n",val2);
            //PAUSE
            val21 = hill->ExpValue(r, GeneIndex, j);		//bechmark value
            add_t = hill->ExpValue(r, GeneIndex, j+1);
            exp_delta = (add_t-val21)/interval_t;
            //printf("cal_exp_delta=%f\n",exp_delta);
            //PAUSE
            val22= hill->CalValue(r, GeneIndex, j);
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
    //PAUSE
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
            
            //i=hill->m_GeneIndex;
            for(j = 0; j < NumTrials; j += time_interval)
            {
                
                cal_profile[r][j]=hill->CalValue(r, GeneIndex, j);
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

