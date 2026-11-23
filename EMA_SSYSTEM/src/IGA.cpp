#include "IGA.h"
#define	Debug_mode  0
extern int test_mode;
IGA::IGA(COrthogonalArray *OAtable,S_System_Divide *ssystem,Knowledge *know,Objfunction *objection)
{
    OA = OAtable;
    s_system = ssystem;
    knowledge = know;
    obj = objection;
    pc = 0;
    ps = 0;
    pm1 =0;
    pm2 = 0;
    NO_GEN =0;
    NO_POP =0;
    NO_g =0;
    NO_h =0;
    PARAM_NUM =0;
    NumGeneParam = 0;
    NumCrossover =0;
    same_interval=false;
    m_size1 = 0;
    m_size2 = 0;
    s_size = 0;
    BestOne = 0;
    NumVars = 0;
    NumRuns = 0;
    NumTrials = 0;
    NumTFs = 0;
    Num_nochange = 0;
    alpha_index = 0;
    beta_index = 0;
    h_start = 0;
    gloc_possible = 0;
    hloc_possible = 0;
    temperature = 0;
    SA = 0;
    cx1 = 0;
    cx2 = 0;
    GeneIndex = 0;
    runtime = 0;
    runnow = 0;
    ptr1 = NULL;
    target = NULL;
    Chromosomes = NULL;
    SortIndex = NULL;
    Lower_bound = NULL;
    Upper_bound = NULL;
    Interval = NULL;
    IntervalChangeByTemplate = NULL;
    sortloc = NULL;
    BestValueEachGeneration = NULL;
    
    NumFitnessFunction = NULL;
    DoMutationIndex = NULL;
    //printf("finish ini iga\n");
}
IGA::~IGA()
{
    int i;
    if(target){delete [] target;}
    if(Chromosomes){delete [] Chromosomes;}
    if(SortIndex){delete [] SortIndex;}
    if(Lower_bound){delete [] Lower_bound;}
    if(Upper_bound){delete [] Upper_bound;}
    if(Interval){delete [] Interval;}
    if(IntervalChangeByTemplate){delete [] IntervalChangeByTemplate;}
    if(sortloc)
    {
        for(i=0;i<2;i++)
            delete [] sortloc[i];
        delete [] sortloc;
    }
    if(BestValueEachGeneration){delete [] BestValueEachGeneration;}
    if(NumFitnessFunction){delete [] NumFitnessFunction;}
    
        
    if(DoMutationIndex){delete [] DoMutationIndex;}
    if(ptr1){fclose(ptr1);}
}
//--------------------------------------------------------------------------------------------
//fitness a is better than fitness b (smaller); NaN is the worst
static inline bool Better(double a,double b)
{
    if(a != a) return false;
    if(b != b) return true;
    return a < b;
}
//a bug somewhere else must not look like a successful run with an empty output
static void iga_fatal(const char *where,int run,int gen)
{
    fprintf(stderr,"IGA internal error (%s): the best fitness got worse, run %d generation %d\n",where,run,gen);
    exit(1);
}
int IGA::LocGene(double code)
{
    int c = (int)code;
    if(c == neg_zero) return 0;
    return (c < 0) ? -c : c;
}
double IGA::LocCode(int gene,bool on)
{
    if(on) return gene;
    return (gene == 0) ? neg_zero : -gene;
}
void IGA::set_parameter(int index,int iteration,int run)
{
    char str1[100];
    pc = DEFAULT_PC;
    ps = DEFAULT_PS;
    pm1 = DEFAULT_PM1;
    pm2 = DEFAULT_PM2;
    if(NO_POP < 2)
    {
        fprintf(stderr,"population size must be at least 2\n");
        exit(1);
    }
    if(NO_GEN < 1)
    {
        fprintf(stderr,"number of generations must be at least 1\n");
        exit(1);
    }
    if(NumVars >= -(neg_zero))
    {
        fprintf(stderr,"too many genes (%d)\n",NumVars);
        exit(1);
    }
    PARAM_NUM = 2*NO_g + 2*NO_h + 3;
    NumGeneParam=2*NumVars + 3;//1 gene(3:alpha + beta + delta_t)
    alpha_index=2*NO_g + 2*NO_h;
    beta_index=2*NO_g + 2*NO_h + 1;
    h_start=2*NO_g;
    //the OA crossover works on pairs of parents: the number of crossovers must be even
    //(it was (int)(0.8*NO_POP), odd for e.g. -P 12, and then read target[NumCrossover])
    NumCrossover = (int)((float)pc*NO_POP);
    NumCrossover -= NumCrossover % 2;
    if(NumCrossover < 2){NumCrossover = 2;}
    m_size1 = (int)(NO_POP * PARAM_NUM * pm1);
    m_size2=(int)(NO_POP * PARAM_NUM * pm2);
    s_size =(int)((float)ps * NO_POP);
    SA=0.997;//退火參數
    GeneIndex = index;
    runtime = run;
    gloc_possible = knowledge->get_gPossible();
    hloc_possible = knowledge->get_hPossible();
    sprintf(str1,"progress%d_step%d.txt",GeneIndex,iteration);
    ptr1=fopen(str1,"w");
    if(ptr1 == NULL)
    {
        fprintf(stderr,"cannot write %s\n",str1);
        exit(1);
    }
    if(Debug_mode)
    {
        printf("pc=%f\n",pc);
        printf("NO_g=%d\n",NO_g);
        printf("possible_g=%d\n",gloc_possible);
        printf("possible_h=%d\n",hloc_possible);
        printf("PARAM_NUM=%d\n",PARAM_NUM);
        printf("GeneIndex=%d\n",GeneIndex);
    }
}
void IGA::initialIGA()
{
    int w;
    //every confirmed regulator needs its own slot
    if(NO_g < knowledge->get_gConfirm() || NO_h < knowledge->get_hConfirm() || NO_g > gloc_possible || NO_h > hloc_possible)
    {
        fprintf(stderr,"inconsistent number of regulators (g: %d slots, %d confirmed, %d possible; h: %d slots, %d confirmed, %d possible)\n",
                NO_g,knowledge->get_gConfirm(),gloc_possible,NO_h,knowledge->get_hConfirm(),hloc_possible);
        exit(1);
    }
    target = new unsigned int[NumCrossover];
    Chromosomes = new double[NO_POP * PARAM_NUM]();
    SortIndex = new int[NO_POP];
    Lower_bound = new double[PARAM_NUM]();
    Upper_bound = new double[PARAM_NUM]();
    Interval = new double[PARAM_NUM]();
    IntervalChangeByTemplate = new double[PARAM_NUM]();
    
    sortloc=new int* [2];
    sortloc[0] = new int[gloc_possible + 1];
    sortloc[1] = new int[hloc_possible + 1];
    for(w=0;w<2;w++)
        memset(sortloc[w],0,(PossibleCount(w)+1)*sizeof(int));
    BestValueEachGeneration = new double[NO_GEN]();//NO_GEN generations are run, see run()
    NumFitnessFunction = new unsigned int[NO_GEN]();
    //printf("run=%d\n",ExecutionRuns);
    //PAUSE
    DoMutationIndex =   new int[PARAM_NUM]();
    OA->Create(PARAM_NUM,Chromosomes,NumVars+1);
    OA->m_Repair = [this](double *c){RepairLocations(c);};
    //OA->DisplayOATable();
}
void IGA::RangeVariable()
{
    int i=0;
    ///==============G_Location & Gij============///
    for(i=0;i<h_start;i++)//
    {
        if(i%2==0)//G_location
        {
            Lower_bound[i]=0;
            Upper_bound[i]=NumTFs-1;//5-1=4
        }
        else//Gij
        {
            Lower_bound[i] = s_system->g_downvalue();
            Upper_bound[i] = s_system->g_upvalue();
            Interval[i]=Upper_bound[i]-Lower_bound[i];
        }
    }
    //===============H_Location & Hij===========//
    for(i=h_start;i<alpha_index;i++)//
    {
        if(i%2==0)//G_location
        {
            Lower_bound[i]=0;
            Upper_bound[i]=NumVars-1;//5-1=4
        }
        else//Hij
        {
            Lower_bound[i] = s_system->h_downvalue();
            Upper_bound[i] = s_system->h_upvalue();
            Interval[i]=Upper_bound[i]-Lower_bound[i];
        }
    }
    //================alpha===================//
    i=alpha_index;
    Lower_bound[i] = s_system->alpha_downvalue();
    Upper_bound[i] = s_system->alpha_upvalue();
    Interval[i]=Upper_bound[i]-Lower_bound[i];
    //================beta====================//
    i=beta_index;
    Lower_bound[i] = s_system->beta_downvalue();
    Upper_bound[i] = s_system->beta_upvalue();
    Interval[i]=Upper_bound[i]-Lower_bound[i];
    //========deltat============///
    i=PARAM_NUM-1;
    //printf("i=%d\n",i);
    Lower_bound[i] = s_system->deltat_downvalue();
    Upper_bound[i] = s_system->deltat_upvalue();
    Interval[i]    = Upper_bound[i]-Lower_bound[i];
    if(Interval[i] == 0){same_interval = true;}
}
//cc of the current best individual (does not change function_cost)
double IGA::BestCC()
{
    unsigned int cost = obj->function_cost;
    obj->fitness(Chromosomes+BestOne*PARAM_NUM);
    obj->function_cost = cost;
    return obj->last_cc;
}
//One execution run: NO_GEN generations (Evaluate, Selection, Crossover, Mutation).
//Early stop: nochange (config file, Num_nochange) is the number of consecutive generations without
//any change of the best fitness after which the run stops; nochange <= 0 switches it off.
//(It used to be `fitness_same <= nochange` together with NO_GEN+1 generations.)
void IGA::run(int run)
{
    
    runnow = run;
    int i;
    int step;
    int NowGeneration;
    int fitness_same;
    double nowfitness=1000000;
    double difference=100;
    double *p;
    int print_no = NO_GEN/10;
    if(print_no < 1){print_no = 1;}//-G < 10
    step = obj->get_time_interval();
    //printf("step:%d\n",step);
    //PAUSE
    BestOne = 0;
    initialpop();
    NowGeneration=0;
    fitness_same=0;
    while((Num_nochange <= 0 || fitness_same < Num_nochange) && (NowGeneration < NO_GEN))
    {
        Evaluate();
        BestValueEachGeneration[NowGeneration] = obj->EVAL[BestOne];
        //printf("Evaluate\tBest=%d\t%e\n",BestOne,obj->EVAL[BestOne]);
        if(NowGeneration>0 && Better(BestValueEachGeneration[NowGeneration-1], obj->EVAL[BestOne]))
        {
            fprintf(stderr,"Bestvalue[%d]=%e\n",NowGeneration-1,BestValueEachGeneration[NowGeneration-1]);
            fprintf(stderr,"EVAL[%d]:%e\n",BestOne,obj->EVAL[BestOne]);
            iga_fatal("Evaluate",runnow,NowGeneration);
        }
        Selection();
        //printf("Selection\tBest=%d\t%e\n",BestOne,obj->EVAL[BestOne]);
        
        if(NowGeneration>0 && Better(BestValueEachGeneration[NowGeneration-1], obj->EVAL[BestOne]))
            iga_fatal("Selection",runnow,NowGeneration);
        Crossover();
        //printf("Crossover\tBest=%d\t%e\n",BestOne,obj->EVAL[BestOne]);
        
        if(NowGeneration>0 && Better(BestValueEachGeneration[NowGeneration-1], obj->EVAL[BestOne]))
            iga_fatal("Crossover",runnow,NowGeneration);
        Mutation();
        //printf("mutation\tBest=%d\t%e\n",BestOne,obj->EVAL[BestOne]);
        if(NowGeneration>0 && Better(BestValueEachGeneration[NowGeneration-1], obj->EVAL[BestOne]))
            iga_fatal("Mutation",runnow,NowGeneration);
        
        ChangeMutationRange();
        NumFitnessFunction[NowGeneration] = obj->function_cost;
        difference = nowfitness-obj->EVAL[BestOne];
        if(difference==0)
            fitness_same++;
        else
            fitness_same=0;
        nowfitness=obj->EVAL[BestOne];
        if(NowGeneration % print_no == 0)
        {
            double cc = BestCC();
            fprintf(ptr1,"%4d %4d %4d %6d %e %1.3f\n",runnow,fitness_same,NowGeneration,obj->function_cost,obj->EVAL[BestOne],cc);
            fflush(ptr1);
            if(test_mode)
            {
                printf("%4d %4d %4d %6d %e %1.3f\n",runnow,fitness_same,NowGeneration,obj->function_cost,obj->EVAL[BestOne],cc);
                p=Chromosomes+BestOne*PARAM_NUM;
                obj->Decode(p);
                obj->printCurrentModelParameters(p,0);
            }
        }
        NowGeneration++;
    }
    p=Chromosomes+BestOne*PARAM_NUM;
    obj->RecordRun(runnow,p);//cc and simulated profile of this run's best model
    obj->Decode(p);//reflash to the best chromosome
    for(i=0;i<NumGeneParam;i++)
        obj->ModelParameters[runnow][i]=obj->GeneParamValues[i];
    obj->ModelParameters[runnow][NumGeneParam]=obj->EVAL[BestOne];
}
//sortloc[w] = the genes that can regulate: the confirmed ('+'/'-') ones first, then the '?' ones
void IGA::BuildSortLoc()
{
    int w,x,k,conf;
    for(w=0;w<2;w++)
    {
        k = 0;
        for(x=0;x<NumVars;x++)
        {
            if(knowledge->get_LocMask(w, x)==0 && knowledge->get_choose(w, x)==1)
            {
                if(k < PossibleCount(w)){sortloc[w][k]=x;}
                k++;
            }
        }
        conf = k;
        for(x=0;x<NumVars;x++)
        {
            if(knowledge->get_LocMask(w, x)==0 && knowledge->get_choose(w, x)==0)
            {
                if(k < PossibleCount(w)){sortloc[w][k]=x;}
                k++;
            }
        }
        if(k != PossibleCount(w) || conf != ConfirmedCount(w))
        {
            fprintf(stderr,"inconsistent knowledge: %d possible regulators found, %d expected\n",k,PossibleCount(w));
            exit(1);
        }
    }
}
//value drawn from the allowed range of the gene (index = chromosome index of Gij/Hij)
double IGA::RandomValue(int w,int gene,int index)
{
    int rg = knowledge->get_range(w, gene);
    if(rg==2)//unknown
        return Lower_bound[index]+Interval[index]*RandomFloat();
    else if(rg==1)//positive
        return 0+Upper_bound[index]*RandomFloat();
    else if(rg==-1)//negative
        return Lower_bound[index]+(0-Lower_bound[index])*RandomFloat();
    return 0;
}
//initial locations of one chromosome (w=0: G, w=1: H)
void IGA::InitLocations(int w,double *chrom)
{
    int no = SlotCount(w), nconf = ConfirmedCount(w), poss = PossibleCount(w), base = SlotBase(w);
    int nfree = no - nconf;
    int s,k,t;
    std::vector<int> pool;
    //the confirmed regulators: dedicated slots, always on
    for(s=0;s<nconf;s++)
        chrom[base+2*s] = LocCode(sortloc[w][s],true);
    if(nfree <= 0)
        return;
    for(k=nconf;k<poss;k++)
        pool.push_back(sortloc[w][k]);
    if(SlotsIdTF(w))
    {
        //one slot for every free gene; a random number (0..nfree) of them is switched off
        int zeros = RandomInt(nfree+1);
        std::vector<int> off(nfree,0);
        for(s=0;s<nfree;s++)
            off[s] = s;
        for(s=0;s<zeros;s++)
        {
            t = s + RandomInt(nfree-s);
            std::swap(off[s],off[t]);
        }
        std::vector<bool> isoff(nfree,false);
        for(s=0;s<zeros;s++)
            isoff[off[s]] = true;
        for(s=0;s<nfree;s++)
            chrom[base+2*(nconf+s)] = LocCode(pool[s],!isoff[s]);
    }
    else
    {
        //random subset: nfree distinct free genes, all on
        for(s=0;s<nfree;s++)
        {
            t = s + RandomInt((int)pool.size()-s);
            std::swap(pool[s],pool[t]);
            chrom[base+2*(nconf+s)] = LocCode(pool[s],true);
        }
    }
}
void IGA::initialpop()
{
    //printf("initial start!!\n");
    int i,j,s,x;
    double *p1;
    int gene;
    obj->function_cost = 0;
    if(Debug_mode)
    {
        for(x=0;x<PARAM_NUM;x++){printf("Low[%d]=%f up[%d]=%f inter[%d]=%f\n",x,Lower_bound[x],x,Upper_bound[x],x,Interval[x]);}
        PAUSE
    }
    BuildSortLoc();
    if(Debug_mode)
    {
        printf("possible G locations:\n");
        for(x=0;x<gloc_possible;x++)
        {
            printf("%d:%d ",x,sortloc[0][x]);
        }
        printf("\npossible H locations:\n");
        for(x=0;x<hloc_possible;x++)
        {
            printf("%d:%d ",x,sortloc[1][x]);
        }
        PAUSE
    }
    
    for(x=0;x<PARAM_NUM;x++){IntervalChangeByTemplate[x]=Interval[x];}
    for(i=0;i<NO_POP;i++)
    {
        p1 = Chromosomes+i*PARAM_NUM;
        InitLocations(0,p1);
        InitLocations(1,p1);
        //Gij / Hij, in the range of the knowledge of their gene
        for(s=0;s<NO_g;s++)
        {
            gene = LocGene(p1[2*s]);
            p1[2*s+1] = knowledge->get_noregulation() ? 0 : RandomValue(0,gene,2*s+1);
        }
        for(s=0;s<NO_h;s++)
        {
            gene = LocGene(p1[h_start+2*s]);
            p1[h_start+2*s+1] = RandomValue(1,gene,h_start+2*s+1);
        }
        for(j=alpha_index;j<PARAM_NUM;j++)//alpha,,beta,deltat
            p1[j]=Lower_bound[j]+Interval[j]*RandomFloat();
    }//end of POPSIZE
    //PAUSE
}
void IGA::ChangeMutationRange()
{
    int i;
    for(i=0;i<PARAM_NUM;i++)
        IntervalChangeByTemplate[i] = Interval[i]* temperature;
    temperature*=SA;
    if(temperature<0.000001)
        temperature=0.01;
}
void IGA::Evaluate()
{
    //printf("Evaluate\n");
    int	i;
    double *p1;
    
    p1=Chromosomes;
    for(i=0;i<NO_POP;i++,p1+=PARAM_NUM)
    {
        obj->EVAL[i] = obj->fitness(p1);
        //printf("%d=%f\n",i,obj->EVAL[i]);
        //PAUSE
    }
    //PAUSE
    
    BestOne = 0;//not the best of the previous generation/run: a NaN there would stay "best"
    for(i=1;i<NO_POP;i++)
    {
        if(Better(obj->EVAL[i], obj->EVAL[BestOne]))
        {
            BestOne = i;
        }
    }
    //printf("Evaluate best index:%d\n",BestOne);
}
void IGA::Selection()
{
    //printf("Selection\n");
    int	i,j,temp;
    
    for(i=0;i<NO_POP;i++)
        SortIndex[i] = i;
    
    // sort by value of fitness function (NaN last)
    for(i=0;i<NO_POP-1;i++)
    {
        for(j=i+1; j<NO_POP;j++)
        {
            if(Better(obj->EVAL[SortIndex[j]], obj->EVAL[SortIndex[i]]))
            {
                temp=SortIndex[j];
                SortIndex[j] = SortIndex[i];
                SortIndex[i] =temp ;
            }
        }
    }
    
    BestOne = SortIndex[0];
    //for(i=0;i<POP_SIZE;i++){printf("%f \n", EVAL[SortIndex[i]]);}
    
    for(i=NO_POP-1;i>=NO_POP-s_size;i--)
    {
        for(j=0;j<PARAM_NUM;j++)
            *(Chromosomes+SortIndex[i]*PARAM_NUM+j) = *(Chromosomes+PARAM_NUM*BestOne+j);
        obj->EVAL[SortIndex[i]]=obj->EVAL[BestOne];
    }
    //printf("Selection best index:%d\n",BestOne);
}

void IGA::Crossover()
{
    //printf("crossover\n");
    int i;
    target[0]=BestOne;//parent1
    for(i=1;i<NumCrossover;i++){
        target[i] = RandomInt(NO_POP);//parent2
    }
    for(i=0;i+1<NumCrossover;i+=2)//NumCrossover is even
    {
        while(target[i]==target[i+1])
        {
            target[i+1] = RandomInt(NO_POP);
        }
        //printf("target[%d]:%f,target[%d]:%f\n",target[i],obj->EVAL[target[i]],target[i+1],obj->EVAL[target[i+1]]);
        CX_mutation(target[i], target[i+1]);
        //printf("crossover: ==========%d %d=========\n",target[i],target[i+1]);
        //PAUSE
        //printf("target[%d]:%f,target[%d]:%f\n",cx1,obj->EVAL[cx1],cx2,obj->EVAL[cx2]);
        //PAUSE
        OA->IntelliegentCrossover(cx1,cx2,BestOne);
        BestOne = OA->set_bestindex();
    }
    //printf("crossover best index:%d\n",BestOne);
}
//New location code of a slot (see IGA.h). Only the free ('?') slots change.
//The new gene is a free gene that is not on in another slot of the chromosome, so the number of
//regulators never drops below the number of slots by creating duplicates.
//Random subset mode (idTF false): the gene is always on (the number of regulators stays NO).
//Mode idTF: the location is switched on/off at random as before.
//A '+'/'-' gene is never touched here, a mutation can neither remove it nor change its sign range.
double IGA::MutateLocation(int w,const double *chrom,int slot)
{
    int nconf = ConfirmedCount(w), poss = PossibleCount(w), base = SlotBase(w), no = SlotCount(w);
    int k,t,gene;
    bool used;
    std::vector<int> cand;
    double oldcode = chrom[base+2*slot];
    if(slot < nconf)
        return LocCode(sortloc[w][slot],true);
    for(k=nconf;k<poss;k++)
    {
        gene = sortloc[w][k];
        used = false;
        for(t=0;t<no && !used;t++)
        {
            if(t != slot && LocOn(chrom[base+2*t]) && LocGene(chrom[base+2*t]) == gene)
                used = true;
        }
        if(!used)
            cand.push_back(gene);
    }
    if(cand.empty())
        return oldcode;
    gene = cand[RandomInt((int)cand.size())];
    if(SlotsIdTF(w))
        return LocCode(gene, (rand() % 2) != 0);
    return LocCode(gene,true);
}
//The crossover combines (location,value) pairs of two parents slot by slot. The dedicated slots of the
//'+'/'-' genes are the same in both parents, but two free slots of the child can hold the same gene
//(then one regulator is lost): move the later one to a gene that is not used.
void IGA::RepairLocations(double *chrom)
{
    int w,s,t,base,nconf,no,gene;
    double code;
    bool dup;
    for(w=0;w<2;w++)
    {
        base = SlotBase(w);
        nconf = ConfirmedCount(w);
        no = SlotCount(w);
        if(w==0 && knowledge->get_noregulation())
            continue;
        for(s=nconf;s<no;s++)
        {
            code = chrom[base+2*s];
            if(!LocOn(code))
                continue;
            gene = LocGene(code);
            dup = false;
            for(t=0;t<s && !dup;t++)
            {
                if(LocOn(chrom[base+2*t]) && LocGene(chrom[base+2*t]) == gene)
                    dup = true;
            }
            if(!dup)
                continue;
            code = MutateLocation(w,chrom,s);
            if(code == chrom[base+2*s])
                continue;//every free gene is used
            chrom[base+2*s] = code;
            if(LocOn(code))
                chrom[base+2*s+1] = RandomValue(w,LocGene(code),base+2*s+1);
        }
    }
}
void IGA::Mutation()
{
    //printf("start mutation\n");
    int i;//nGene;
    int m_pop,m_param,point;
    int g_index;
    int h_index;
    int gij_index;
    int hij_index;
    int slot,nfree;
    double gij,hij;
    double value, oldvalue, newvalue;
    double *chrom;
    
    if(!knowledge->get_noregulation() && NO_g > 0)
    {
        //============g_location====================================//
        nfree = NO_g - knowledge->get_gConfirm();
        if(nfree > 0 && !knowledge->get_fixg())
        {
            for(i=0;i<m_size2;i++)
            {
                do{
                    m_pop = RandomInt(NO_POP);
                }while(m_pop == BestOne);
                chrom = Chromosomes+PARAM_NUM*m_pop;
                slot = knowledge->get_gConfirm() + RandomInt(nfree);
                g_index = 2*slot;
                oldvalue = chrom[g_index];
                newvalue = MutateLocation(0,chrom,slot);
                chrom[g_index] = newvalue;
                //every time the location changes the Gij has to follow (new gene = new sign range)
                if(newvalue != oldvalue)
                {
                    gij_index = g_index + 1;
                    gij = chrom[gij_index];
                    chrom[gij_index] = MutateOneValue(gij, gij_index, static_cast<int>(newvalue));
                }
            }
        }
        //=============Gij========================================//
        for(i=0;i<m_size1;i++)
        {
            do{
                m_pop = RandomInt(NO_POP);
            }while(m_pop == BestOne);
            gij_index = 2*RandomInt(NO_g)+1;
            g_index = static_cast<int>(*(Chromosomes+PARAM_NUM*m_pop+(gij_index-1)));
            gij = *(Chromosomes+PARAM_NUM*m_pop+gij_index);
            *(Chromosomes+PARAM_NUM*m_pop+gij_index)=MutateOneValue(gij,gij_index,g_index);
        }
        //printf("finish Gij\n");
    }
    
    
    
    //============h_location=================================//
    if(NO_h > 0)//no Hij is searched when all Hij are fixed to 0
    {
        nfree = NO_h - knowledge->get_hConfirm();
        if(nfree > 0 && !knowledge->get_fixh())
        {
            for(i=0;i<m_size2;i++)
            {
                do{
                    m_pop = RandomInt(NO_POP);
                }while(m_pop == BestOne);
                chrom = Chromosomes+PARAM_NUM*m_pop;
                slot = knowledge->get_hConfirm() + RandomInt(nfree);
                h_index = h_start + 2*slot;
                oldvalue = chrom[h_index];
                newvalue = MutateLocation(1,chrom,slot);
                chrom[h_index] = newvalue;
                if(newvalue != oldvalue)
                {
                    hij_index = h_index + 1;
                    hij = chrom[hij_index];
                    chrom[hij_index] = MutateOneValue(hij, hij_index, static_cast<int>(newvalue));
                }
            }
        }
        //============hij=================================//
        for(i=0;i<m_size1;i++)
        {
            do{
                m_pop = RandomInt(NO_POP);
            }while(m_pop == BestOne);
            hij_index = h_start + 2*RandomInt(NO_h)+1;
            h_index = static_cast<int>(*(Chromosomes+PARAM_NUM*m_pop+(hij_index-1)));
            hij = *(Chromosomes+PARAM_NUM*m_pop+hij_index);
            *(Chromosomes+PARAM_NUM*m_pop+hij_index)=MutateOneValue(hij,hij_index,h_index);
        }
    }//end of NO_h > 0
    //============alpha, beta==============================================================
    for(i=0;i<m_size1;i++)
    {
        do{
            m_pop = RandomInt(NO_POP);
        }while(m_pop == BestOne);//choose the chromosome to mutation
        point=rand() % 2;
        m_param = alpha_index + point;
        value=*(Chromosomes+PARAM_NUM*m_pop+m_param);
        //printf("basel transcription and degrade value=%f\n",value);
        *(Chromosomes+PARAM_NUM*m_pop+m_param)=MutateOneValue(value,m_param,0);
    }
    //printf("finish others\n");
    if(!same_interval)
    {
        for(i=0;i<m_size1;i++)//delta_t
        {
            do{
                m_pop = RandomInt(NO_POP);
            }while(m_pop == BestOne);//choose the chromosome to mutation
            m_param = PARAM_NUM-1;
            value=*(Chromosomes+PARAM_NUM*m_pop+m_param);
            *(Chromosomes+PARAM_NUM*m_pop+m_param)=MutateOneValue(value,m_param,0);
        }
    }
    //PAUSE
    //printf("mutation best index:%d\n",BestOne);
    //printf("mutant finish\n");
}
double IGA::MutateOneValue(double oldValue,int index,int fixloc)
{
    
    //printf("oldValue=%f\n",oldValue);
    //printf("index=%d\n",index);
    //printf("fixloc=%d\n",fixloc);
        //PAUSE
    double value,r, gradient;
    //printf("index=%d\n",index);
    //PAUSE
    value = 0;
    if(index < h_start)
    {
        if(index % 2 == 0)//g_Location: mutated by MutateLocation, never here
        {
            value = oldValue;
        }
        else//Gij
        {
            
            //printf("mutate gij\n");
            //printf("before mutant:%f, range[%d]=%d\n",oldValue,fixloc,knowledge->get_range(0, fixloc));
            //PAUSE
            if(fixloc < 0)
                value = oldValue;
            else
            {
                if(knowledge->get_range(0, fixloc) == 2)
                {
                    //printf("? mutant old value:%f\n",oldValue);
                    //PAUSE
                    do
                    {
                        r = RandomFloat();
                    }while(r == 0.5);
                    gradient = 0.5 * IntervalChangeByTemplate[index] * tan(PI*r);
                    value = oldValue + gradient;
                    if((value > Upper_bound[index])||(value < Lower_bound[index]))
                    {
                        value = Lower_bound[index] + Interval[index] * RandomFloat();
                    }
                    
                }
                else if(knowledge->get_range(0, fixloc) == 1)
                {
                    //printf("positive\n");
                    //PAUSE
                    do
                    {
                        r = RandomFloat();
                    }while(r == 0.5);
                    gradient = 0.5 * IntervalChangeByTemplate[index] * tan(PI*r);
                    value = oldValue + gradient;
                    if((value > Upper_bound[index])||(value < 0))
                    {
                        value = 0 + Upper_bound[index] * RandomFloat();
                    }
                }
                else if(knowledge->get_range(0, fixloc)==-1)
                {
                    //printf("negative\n");
                    //PAUSE
                    do
                    {
                        r = RandomFloat();
                    }while(r == 0.5);
                    gradient = 0.5 * IntervalChangeByTemplate[index] * tan(PI*r);
                    value = oldValue + gradient;
                    if((value > 0)||(value < Lower_bound[index]))
                    {
                        value= Lower_bound[index] + (0-Lower_bound[index]) * RandomFloat();
                    }
                    
                }
            }
            
            //printf("after mutant beta:%f\n",value);
            //PAUSE
        }
        
    }
    else if(index>=h_start && index<alpha_index)
    {
        if(index % 2 == 0)//h_Location: mutated by MutateLocation, never here
        {
            value = oldValue;
        }
        else//Hij
        {
            
            //printf("mutate hij\n");
            //printf("before mutant:%f, range[%d]=%d\n",oldValue,fix,range[fix]);
            //PAUSE
            if(fixloc < 0)
                value = oldValue;
            else
            {
                if(knowledge->get_range(1, fixloc) == 2)
                {
                    //printf("? mutant old value:%f\n",oldValue);
                    //PAUSE
                    do
                    {
                        r = RandomFloat();
                    }while(r == 0.5);
                    gradient = 0.5 * IntervalChangeByTemplate[index] * tan(PI*r);
                    value = oldValue + gradient;
                    if((value > Upper_bound[index])||(value < Lower_bound[index]))
                    {
                        value = Lower_bound[index] + Interval[index] * RandomFloat();
                    }
                    
                    
                }
                else if(knowledge->get_range(1, fixloc) == 1)
                {
                    //printf("positive\n");
                    //PAUSE
                    do
                    {
                        r = RandomFloat();
                    }while(r == 0.5);
                    gradient = 0.5 * IntervalChangeByTemplate[index] * tan(PI*r);
                    value = oldValue + gradient;
                    if((value > Upper_bound[index])||(value < 0))
                    {
                        value = 0 + Upper_bound[index] * RandomFloat();
                    }
                }
                else if(knowledge->get_range(1, fixloc) == -1)
                {
                    //printf("negative\n");
                    //PAUSE
                    do
                    {
                        r = RandomFloat();
                    }while(r == 0.5);
                    gradient = 0.5 * IntervalChangeByTemplate[index] * tan(PI*r);
                    value = oldValue + gradient;
                    if((value > 0)||(value < Lower_bound[index]))
                    {
                        value= Lower_bound[index] + (0-Lower_bound[index]) * RandomFloat();
                    }
                    
                }
            }
            
            //printf("after mutant beta:%f\n",value);
            //PAUSE
        }
        
    }
    else if (index == PARAM_NUM-1)    //deltat
    {
        //printf("delta_t mutant\n");
        if(same_interval){value = oldValue;}
        else
        {
            do
            {
                r = RandomFloat();
            }while(r == 0.5);
            gradient = 0.5 * IntervalChangeByTemplate[index] * tan(PI*r);
            value = oldValue + gradient;
            if((value > Upper_bound[index])||(value < Lower_bound[index]))
            {
                value = Lower_bound[index] + Interval[index] * RandomFloat();
            }

        }
    }
    
    else	//alpha,beta
    {
        //printf("mutant alpha,beta\n");
        do
        {
            r = RandomFloat();
        }while(r == 0.5);
        gradient = 0.5 * IntervalChangeByTemplate[index] * tan(PI*r);
        value = oldValue + gradient;
        if((value > Upper_bound[index])||(value < Lower_bound[index]))
        {
            value = Lower_bound[index] + Interval[index] * RandomFloat();
        }

    }
    //printf("after mutate value=%f\n",value);
    //PAUSE
    return value;
}
void IGA::CX_mutation(unsigned int index_one,unsigned int index_two)
{
    unsigned int i;
    
    //double value;
    bool reCom=false;
    int Index[2];
    //int DoMutationIndex[PARAM_NUM];
    int loc;
    int w, slot;
    double oldcode, newcode;
    if(!(Better(obj->EVAL[index_two], obj->EVAL[index_one])))//index_one is the better (or equal) one
    {
        //printf("EVAL[%d]=%f EVAL[%d]=%f\n",*index_one,EVAL[*index_one],*index_two,EVAL[*index_two]);
        Index[0]=index_one;
        Index[1]=index_two;
    }
    else
    {
        Index[0]=index_two;
        Index[1]=index_one;
    }
    cx1=Index[0];
    cx2=Index[1];
    //printf("Index[0]=%d, Index[1]=%d\n",Index[0],Index[1]);
    //printf("%d,%d\n",index_one,index_two);
    //PAUSE
    OA->m_Parents[0]=Chromosomes+Index[0]*PARAM_NUM;
    OA->m_Parents[1]=Chromosomes+Index[1]*PARAM_NUM;
    if(Debug_mode)
    {
        printf("\nIndex0:%d fitness:%f\n",Index[0],obj->EVAL[Index[0]]);
        for(i=0;i<PARAM_NUM;i++)
        {
            printf("%f ",OA->m_Parents[0][i]);
        }
        printf("\n");
        printf("Index1:%d fitness:%f\n",Index[1],obj->EVAL[Index[1]]);
        for(i=0;i<PARAM_NUM;i++)
        {
            printf("%f ",OA->m_Parents[1][i]);
        }
        printf("\n");
        PAUSE
    }
    OA->m_OATable[0][PARAM_NUM]=obj->EVAL[Index[0]];
    OA->m_OATable[OA->m_Height][PARAM_NUM]=obj->EVAL[Index[1]];
    for(i=0;i<PARAM_NUM;i++)
        DoMutationIndex[i]=0;
    if(same_interval){DoMutationIndex[PARAM_NUM-1]=1;}
    for(i=0;i<PARAM_NUM;i++)
    {
        OA->m_OATable[0][i]=OA->m_Parents[0][i];
        OA->m_OATable[OA->m_Height][i]=OA->m_Parents[1][i];
        if(OA->m_OATable[0][i]==OA->m_OATable[OA->m_Height][i] ) //&& i < m_HalfNumParameters
        {
            //printf("crossover the same\n");
            //PAUSE
            if(DoMutationIndex[i]==0 )
            {
                
                if(i < (unsigned int)alpha_index)
                {
                    w = ((int)i < h_start) ? 0 : 1;
                    if(w == 0 && knowledge->get_noregulation())
                        continue;
                    if((int)i % 2 == 0)//g_Location / h_Location
                    {
                        slot = ((int)i - SlotBase(w))/2;
                        oldcode = OA->m_Parents[1][i];
                        newcode = MutateLocation(w,OA->m_Parents[1],slot);
                        OA->m_Parents[1][i] = newcode;
                        OA->m_OATable[OA->m_Height][i]=OA->m_Parents[1][i];
                        DoMutationIndex[i]=1;
                        //the value of the Gij/Hij has to follow EVERY change of the location
                        //(the old flag was never reset, so after the first unchanged location it stopped doing so)
                        if(!SlotsFixed(w) && newcode != oldcode)
                        {
                            OA->m_Parents[1][i+1] = MutateOneValue(OA->m_Parents[1][i+1],i+1,static_cast<int>(newcode));
                            OA->m_OATable[OA->m_Height][i+1]=OA->m_Parents[1][i+1];
                            DoMutationIndex[i+1]=1;
                        }
                    }
                    else//Gij / Hij
                    {
                        if(DoMutationIndex[i] == 0)
                        {
                            loc = static_cast<int>(OA->m_Parents[1][i-1]);
                            OA->m_Parents[1][i]=MutateOneValue(OA->m_Parents[1][i],i,loc);
                            OA->m_OATable[OA->m_Height][i]=OA->m_Parents[1][i];
                            DoMutationIndex[i]=1;
                        }
                    }
                    reCom=true;
                }
                else
                {//mutate b,n,k,degrade
                    //printf("mutate b,n,k,degrade\n");
                    OA->m_Parents[1][i]=MutateOneValue(OA->m_Parents[1][i],i,0);
                    OA->m_OATable[OA->m_Height][i]=OA->m_Parents[1][i];
                    DoMutationIndex[i]=1;
                    reCom=true;
                }
            }
        }
    }
    if(reCom)
    {
        obj->EVAL[Index[1]]=obj->fitness(OA->m_Parents[1]);
        //double testvalue = fitness(m_OATable[m_Height]);
        //printf("%f, %f\n",EVAL[Index[1]],testvalue);
        //PAUSE
        /*
         printf("after mutant parant1 fitness:%f\n",EVAL[Index[1]]);
         for(i=0;i<m_NumParameters;i++)
         {
         printf("%f ",m_Parents[1][i]);
         }
         printf("\n");
         PAUSE
         */
        //EVAL[Index[1]]=fitness(m_Parents[1]);//caculate the new fitness after mutation
        OA->m_OATable[OA->m_Height][PARAM_NUM]=obj->EVAL[Index[1]];
    }
    //printf("domutation over!!\n");
    //PAUSE
}
