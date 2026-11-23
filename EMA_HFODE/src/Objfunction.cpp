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
    last_cc = 0;
    fitness_type = 0;
    GeneIndex = 0;
    cc_mode = 0;
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
    memset(GeneParamValues, 0, NumGeneParam*sizeof(double));
    cal_profile = new double* [NumRuns];
    for(i=0;i<NumRuns;i++)
        cal_profile[i] = new double[NumTrials]();//value-initialised (zero)
    NumConnections = No_reg;
    time_interval = (time_step < 1) ? 1 : time_step;//a step of 0 would never advance the time loops
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
//Decode one chromosome row into GeneParamValues
//(b, TransMax, n(0..NumVars-1), k(0..NumVars-1), degrade).
//Chromosome layout: loc(s)|n(s) for slot s=0..C-1, k(s) for s=0..C-1, b, TransMax, degrade, mask.
//Slot s is active when bit s of the mask is set; a negative location -g means "gene g is
//switched off" (n and k of g are 0).  Confirmed regulators are the first slots with the
//mask bit always set, so they are never masked out or switched off.
void Objfunction::DecodeRow(const double *row)
{
    int i, pp, s;
    for(i=0;i<NumGeneParam;i++){
        GeneParamValues[i] = 0;//initial
    }
    GeneParamValues[0] = row[PARAM_NUM-4];//bi
    GeneParamValues[1] = row[PARAM_NUM-3];//TransMax
    GeneParamValues[NumGeneParam-1] = row[PARAM_NUM-2];//degrade
    unsigned long long MM = (unsigned long long)row[PARAM_NUM-1];//Mask (64 bit, up to MAX_SLOTS slots)
    for(s=0;s<NumConnections;s++)
    {
        pp = (int)row[2*s];
        if(pp >= 0)
        {
            if(pp >= NumTFs) continue;
            if((MM >> s) & 1ULL)
            {
                GeneParamValues[2+pp] = row[2*s+1];//n
                GeneParamValues[2+NumVars+pp] = row[2*NumConnections+s];//k
            }
        }else{
            if(pp == neg_zero)
                pp = 0;
            else
                pp = abs(pp);
            if(pp >= NumTFs) continue;
            GeneParamValues[2+pp] = 0;
            GeneParamValues[2+NumVars+pp] = 0;
        }
    }
}
void Objfunction::Decode(int PopIndex)
{
    DecodeRow(Chromosomes[PopIndex]);
}
double Objfunction::OAfitenss(double **OATable, int OATableIndex)
{
    DecodeRow(OATable[OATableIndex]);
    return fitnessValue(GeneParamValues);
}
double Objfunction::fitness(int PopIndex)
{
    Decode(PopIndex);
    return fitnessValue(GeneParamValues);
}
double Objfunction::EvaluateBest(int PopIndex)
{
    int r, j;
    double value;
    Decode(PopIndex);
    value = fitnessValue(GeneParamValues, false);
    for(r = 0; r < NumRuns; r++)
        for(j = 0; j < NumTrials; j += time_interval)
            cal_profile[r][j] = hill->CalValue(r, GeneIndex, j);
    best_cc = last_cc;
    return value;
}
double Objfunction::fitnessValue(double *OneChrom, bool count_cost)//(GeneParamValues)
{
    int j, r, x, npts;
    double fitness, val,val2,val21,val22;
    double exp_avg,cal_avg,exp_total,cal_total,cal_diff,exp_diff,cal_diff_pow,exp_diff_pow,cal_exp_multi,c_c,c_c_fitness,add_t,exp_delta,interval_t;
    double sum, totalSum, ErrorValue, denom, cc_sum;
    double temp;
    double max_cc = 1;
    //b
    hill->SetCalb(GeneIndex, OneChrom[0]);
    //TransMax
    hill->SetCalTransMax(GeneIndex, OneChrom[1]);
    //Nij
    for(j = 0; j < NumTFs; j++)
        hill->SetCaln(GeneIndex, j, OneChrom[j+2]);
    //Kij
    for(j = 0; j < NumTFs; j++)
        hill->SetCalk(GeneIndex, j, OneChrom[2+NumVars+j]);
    //degrade[i]
    hill->SetCaldegrade(GeneIndex, OneChrom[NumGeneParam-1]);
    interval_t=hill->GetDeltaT();
    
    hill->CalAllCalValues(GeneIndex); //calculate the estimate profile
    
    //Expression Error
    totalSum = 0;
    cc_sum = 0;
    for(r = 0; r < NumRuns; r++)
    {
        sum = 0;
        exp_total=0;
        cal_total=0;
        npts = 0;//number of time points that are summed (every time_interval-th one)
        for(x = 0; x < NumTrials; x += time_interval)
        {
            exp_total += hill->ExpValue(r, GeneIndex, x);
            cal_total += hill->CalValue(r, GeneIndex, x);
            npts++;
        }
        //the means must be divided by the number of summed points (not NumTrials)
        exp_avg = exp_total/npts;
        cal_avg = cal_total/npts;
        cal_exp_multi = 0;
        exp_diff_pow = 0;
        cal_diff_pow = 0;
        for(j = 0; j < NumTrials; j += time_interval)
        {
            val21 = hill->ExpValue(r, GeneIndex, j);		//bechmark value
            val22 = hill->CalValue(r, GeneIndex, j);
            //--------cc--------------------------//
            exp_diff = val21 - exp_avg;
            cal_diff = val22 - cal_avg;
            cal_exp_multi += exp_diff*cal_diff;
            exp_diff_pow += pow(exp_diff,2);
            cal_diff_pow += pow(cal_diff,2);
            //--------cc---------------------------//
            //The derivative term needs the next observed point.  At the last time point
            //there is none (X[NumTrial] is never observed), so the derivative term is
            //skipped there while the value-error term is still used.
            bool has_next = (j + 1 < NumTrials);
            temp = 0;
            if(has_next && fitness_type >= 2)
            {
                val2 = hill->CalDeltaValue(r, GeneIndex, j);
                add_t = hill->ExpValue(r, GeneIndex, j+1);
                exp_delta = (add_t-val21)/interval_t;
                temp = fabs(exp_delta - val2);
            }
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
        //correlation coefficient of this replicate; undefined (0/0) for a constant
        //profile, then cc = 0 (so the cc fitness term is 1)
        denom = sqrt(cal_diff_pow)*sqrt(exp_diff_pow);
        if(denom > 1e-12)
            c_c = cal_exp_multi/denom;
        else
            c_c = 0;
        if(c_c > 1){c_c = 1;}
        if(c_c < -1){c_c = -1;}
        if(c_c != c_c){c_c = 0;}
        cc_sum += c_c;
    }//end of NumRuns
    
    c_c = cc_sum/(double)NumRuns;//mean cc over all replicates
    c_c_fitness = max_cc - c_c;
    totalSum /= (double)NumRuns;
    if(cc_mode)
        fitness = totalSum + c_c_fitness;
    else
        fitness = totalSum;
    //keep every fitness finite (NaN/inf -> worst) so all GA comparisons are well defined
    if(!(fitness < WORST_FITNESS))
        fitness = WORST_FITNESS;
    last_cc = c_c;
    if(count_cost)
        function_cost ++;
    return fitness;
}
