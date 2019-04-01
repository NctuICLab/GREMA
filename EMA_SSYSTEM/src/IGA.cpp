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
    temperature = 0;
    SA = 0;
    cx1 = 0;
    cx2 = 0;
    GeneIndex = 0;
    runnow = 0;
    target = NULL;
    Chromosomes = NULL;
    SortIndex = NULL;
    Lower_bound = NULL;
    Upper_bound = NULL;
    Interval = NULL;
    IntervalChangeByTemplate = NULL;
    Loc_tmp = NULL;
    eachpop_loc = NULL;
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
    if(Loc_tmp)
    {
        for(i=0;i<2;i++)
            delete [] Loc_tmp[i];
        delete [] Loc_tmp;
    }
    if(eachpop_loc)
    {
        for(i=0;i<2;i++)
            delete [] eachpop_loc[i];
        delete [] eachpop_loc;
    }
    if(sortloc)
    {
        for(i=0;i<2;i++)
            delete [] sortloc[i];
        delete [] sortloc;
    }
    if(BestValueEachGeneration){delete [] BestValueEachGeneration;}
    if(NumFitnessFunction){delete [] NumFitnessFunction;}
    
        
    if(DoMutationIndex){delete [] DoMutationIndex;}
    fclose(ptr1);
}
void IGA::set_parameter(int index,int iteration,int run)
{
    char str1[100];
    pc = DEFAULT_PC;
    ps = DEFAULT_PS;
    pm1 = DEFAULT_PM1;
    pm2 = DEFAULT_PM2;
    PARAM_NUM = 2*NO_g + 2*NO_h + 3;
    NumGeneParam=2*NumVars + 3;//1 gene(3:alpha + beta + delta_t)
    alpha_index=2*NO_g + 2*NO_h;
    beta_index=2*NO_g + 2*NO_h + 1;
    h_start=2*NO_g;
    NumCrossover = (int)((float)pc*NO_POP);
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
    int i;
    target = new unsigned int[NumCrossover];
    Chromosomes = new double[NO_POP * PARAM_NUM];
    SortIndex = new int[NO_POP];
    Lower_bound = new double[PARAM_NUM];
    Upper_bound = new double[PARAM_NUM];
    Interval = new double[PARAM_NUM];
    IntervalChangeByTemplate = new double[PARAM_NUM];
    
    Loc_tmp=new int* [2];//initial selected location
    eachpop_loc=new int* [2];
    for(i=0;i<2;i++)
    {
        Loc_tmp[i] = new int[NumVars];
        eachpop_loc[i] = new int[NumVars];
    }
    sortloc=new int* [2];
    sortloc[0] = new int[gloc_possible];
    sortloc[1] = new int[hloc_possible];
    BestValueEachGeneration = new double[NO_GEN];
    memset(BestValueEachGeneration,0,NO_GEN);
    NumFitnessFunction = new unsigned int[NO_GEN];
    memset(NumFitnessFunction,0,NO_GEN);
    //printf("run=%d\n",ExecutionRuns);
    //PAUSE
    DoMutationIndex =   new int[PARAM_NUM];
    OA->Create(PARAM_NUM,Chromosomes,NumVars+1);
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
    step = obj->get_time_interval();
    //printf("step:%d\n",step);
    //PAUSE
    initialpop();
    NowGeneration=0;
    fitness_same=0;
    while((fitness_same <= Num_nochange) && (NowGeneration <= NO_GEN))
    {
        Evaluate();
        BestValueEachGeneration[NowGeneration] = obj->EVAL[BestOne];
        //printf("Evaluate\tBest=%d\t%e\n",BestOne,obj->EVAL[BestOne]);
        if(NowGeneration>0 && BestValueEachGeneration[NowGeneration-1]< obj->EVAL[BestOne])
        {
            printf("Bestvalue[%d]=%e\n",NowGeneration-1,BestValueEachGeneration[NowGeneration-1]);
            printf("EVAL[%d]:%e\n",BestOne,obj->EVAL[BestOne]);
            fprintf(stderr,"Evaluate error...........");
            exit(0);
        }
        Selection();
        //printf("Selection\tBest=%d\t%e\n",BestOne,obj->EVAL[BestOne]);
        
        if(NowGeneration>0 && BestValueEachGeneration[NowGeneration-1]<obj->EVAL[BestOne])
        {
            fprintf(stderr,"Selection error...........");
            exit(0);
        }
        Crossover();
        //printf("Crossover\tBest=%d\t%e\n",BestOne,obj->EVAL[BestOne]);
        
        if(NowGeneration>0 && BestValueEachGeneration[NowGeneration-1]<obj->EVAL[BestOne])
        {
            fprintf(stderr,"Crossover error...........");
            exit(0);
        }
        Mutation();
        //printf("mutation\tBest=%d\t%e\n",BestOne,obj->EVAL[BestOne]);
        if(NowGeneration>0 && BestValueEachGeneration[NowGeneration-1]<obj->EVAL[BestOne])
        {
            fprintf(stderr,"Mutation error...........");
            exit(0);
        }
        
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
            fprintf(ptr1,"%4d %4d %4d %6d %e %1.3f\n",runnow,fitness_same,NowGeneration,obj->function_cost,obj->EVAL[BestOne],obj->best_cc);
            fflush(ptr1);
            if(test_mode)
            {
                printf("%4d %4d %4d %6d %e %1.3f\n",runnow,fitness_same,NowGeneration,obj->function_cost,obj->EVAL[BestOne],obj->best_cc);
                p=Chromosomes+BestOne*PARAM_NUM;
                obj->Decode(p);
                obj->printCurrentModelParameters(p,0);
            }
        }
        if(NowGeneration>0 && BestValueEachGeneration[NowGeneration-1]<obj->EVAL[BestOne])
        {
            fprintf(stderr,"IGA error...........");
            exit(0);
        }
        NowGeneration++;
    }
    p=Chromosomes+BestOne*PARAM_NUM;
    obj->Decode(p);//reflash to the best chromosome
    for(i=0;i<NumGeneParam;i++)
        obj->ModelParameters[runnow][i]=obj->GeneParamValues[i];
    obj->ModelParameters[runnow][NumGeneParam]=obj->EVAL[BestOne];
}
void IGA::initialpop()
{
    //printf("initial start!!\n");
    int i,j,x;
    int g,h;
    double *p1;
    double value;
    int g_possible_zero;
    int h_possible_zero;
    int sel_gloc,sel_hloc;
    int gloc_index,hloc_index;
    int tmp_loc;
    temperature = 1.0;
    obj->function_cost = 0;
    obj->best_cc=-1;
    value = 0;
    g = 0;
    h = 0;
    if(Debug_mode)
    {
        for(x=0;x<PARAM_NUM;x++){printf("Low[%d]=%f up[%d]=%f inter[%d]=%f\n",x,Lower_bound[x],x,Upper_bound[x],x,Interval[x]);}
        PAUSE
    }
    for(x=0;x<NumVars;x++)
    {
        if(knowledge->get_LocMask(0, x)==0)
        {
            sortloc[0][g]=x;
            g++;
        }
    }
    for(x=0;x<NumVars;x++)
    {
        if(knowledge->get_LocMask(1, x)==0)
        {
            sortloc[1][h]=x;
            h++;
        }
    }
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
    p1=Chromosomes;
    for(i=0;i<NO_POP;i++)
    {
        sel_gloc=0;
        sel_hloc=0;
        gloc_index=0;
        hloc_index=0;
        for(x=0;x<NumTFs;x++){Loc_tmp[0][x]=knowledge->get_LocMask(0, x);}
        for(x=0;x<NumTFs;x++){eachpop_loc[0][x]=knowledge->get_LocMask(0, x);}
        for(x=0;x<NumVars;x++){Loc_tmp[1][x]=knowledge->get_LocMask(1, x);}
        for(x=0;x<NumVars;x++){eachpop_loc[1][x]=knowledge->get_LocMask(1, x);}
        if(!(knowledge->get_fixg()))
        {
            if(Debug_mode)
            {
             printf("Gij original:\n");
             for(x=0;x<NumVars;x++)
                 printf("%d:%d ",x,eachpop_loc[0][x]);
            }
            g_possible_zero = RandomInt(NO_g - knowledge->get_gConfirm()+1);
            if(Debug_mode)
            {
                printf("\n Gij number of zero:%d\n",g_possible_zero);
                PAUSE
            }
            for(x=0;x<g_possible_zero;x++)
            {
                tmp_loc=0;
                do{
                    tmp_loc = RandomInt(NumTFs);
                }while(eachpop_loc[0][tmp_loc]==1||knowledge->get_choose(0, tmp_loc)==1);
                eachpop_loc[0][tmp_loc]=1;
                //printf("%d:mutant loc:%d\n",k,tmp_loc);
            }
            
            if(Debug_mode)
            {
             printf("after mutate Gij:\n");
             for(x=0;x<NumVars;x++)
                 printf("%d:%d ",x,eachpop_loc[0][x]);
             printf("\n");
            }
        }
        if(!(knowledge->get_fixh()))
        {
            if(Debug_mode)
            {
             printf("\nHij original:\n");
             for(x=0;x<NumVars;x++)
                 printf("%d:%d ",x,eachpop_loc[1][x]);
            }
            h_possible_zero = RandomInt(NO_h - knowledge->get_hConfirm()+1);
            if(Debug_mode)
            {
                printf("\n Hij number of zero:%d\n",h_possible_zero);
                PAUSE
            }
            for(x=0;x<h_possible_zero;x++)
            {
                tmp_loc=0;
                do{
                    tmp_loc = RandomInt(NumVars);
                }while(eachpop_loc[1][tmp_loc]==1||knowledge->get_choose(1, tmp_loc)==1);
                eachpop_loc[1][tmp_loc]=1;
                //printf("%d:mutant loc:%d\n",k,tmp_loc);
            }
            if(Debug_mode)
            {
             printf("\nafter mutate Hij:\n");
             for(x=0;x<NumVars;x++)
                 printf("%d:%d ",x,eachpop_loc[1][x]);
             printf("\n");
            }
        }
        
        //printf("%d:",i);
        for(j=0;j<PARAM_NUM;j++,p1++)
        {
            
            if(j<h_start) //g_location and Gij
            {
                if(j%2 == 0)//g_location
                {
                    if(knowledge->get_fixg())
                    {
                        if(knowledge->get_noregulation())
                        {
                            value = 0;
                        }
                        else
                        {
                            sel_gloc = sortloc[0][gloc_index];
                            value = sortloc[0][gloc_index];
                        }
                        
                    }
                    else
                    {
                        if(knowledge->get_idgTF())
                        {
                            tmp_loc = sortloc[0][gloc_index];
                            //printf("\n%d %d\n",gloc_index,tmp_loc);
                            //PAUSE
                            
                            sel_gloc = tmp_loc;//select location
                            if(eachpop_loc[0][sel_gloc] == 1)//mutant to negative
                            {
                                if(tmp_loc == 0)
                                    value = neg_zero;
                                else
                                    value = -(tmp_loc);
                            }
                            else
                                value = tmp_loc;
                        
                        }
                        else
                        {
                            do{
                                tmp_loc = RandomInt(NumVars);
                            }while(Loc_tmp[0][tmp_loc]==1);
                            Loc_tmp[0][tmp_loc]=1;
                            sel_gloc = tmp_loc;
                            value = tmp_loc;
                        }
                        
                    }
                    
                    //printf("g_Location: %f\n",value);
                    //PAUSE
                    gloc_index++;
                    //printf("location initial,knowledge=%d\n",g_NumConnections);
                }
                else//Gij
                {
                    //printf("sel_gloc=%d, j=%d\n",sel_gloc,j);
                    if(knowledge->get_noregulation())
                        value = 0;
                    else
                    {
                        if(knowledge->get_range(0, sel_gloc)==2)
                        {
                            //printf("unknown ");
                            //PAUSE
                            value=Lower_bound[j]+Interval[j]*RandomFloat();
                        }
                        else if(knowledge->get_range(0, sel_gloc)==1)
                        {
                            //printf("positive ");
                            //PAUSE
                            value=0+Upper_bound[j]*RandomFloat();
                        }
                        else if(knowledge->get_range(0, sel_gloc)==-1)
                        {
                            //printf("negative ");
                            //PAUSE
                            value=Lower_bound[j]+(0-Lower_bound[j])*RandomFloat();
                        }
                        //printf("Gij: %f\n",value);
                    }
                }
            }
            else if(j>=h_start && j<alpha_index) //h_location and Hij
            {
                if(j%2 == 0)//h_location
                {
                    if(knowledge->get_fixh())
                    {
                        sel_hloc = sortloc[1][hloc_index];
                        value = sortloc[1][hloc_index];
                    }
                    else
                    {
                        if(knowledge->get_idhTF())
                        {
                            tmp_loc = sortloc[1][hloc_index];
                            //printf("\n%d %d\n",hloc_index,tmp_loc);
                            //PAUSE
                            sel_hloc = tmp_loc;//select location
                            if(eachpop_loc[1][sel_hloc] == 1)//mutant to negative
                            {
                                if(tmp_loc == 0)
                                    value = neg_zero;
                                else
                                    value = -(tmp_loc);
                            }
                            else
                                value = tmp_loc;
                        
                        }
                        else
                        {
                            do{
                                tmp_loc = RandomInt(NumVars);
                            }while(Loc_tmp[1][tmp_loc]==1);
                            Loc_tmp[1][tmp_loc]=1;
                            sel_hloc = tmp_loc;
                            value = tmp_loc;
                        }
                        
                    }
                    //printf("h_Location: %f\n",value);
                    //PAUSE
                    hloc_index++;
                    //printf("Location: %f\n",value);
                    //printf("location initial,knowledge=%d\n",g_NumConnections);
                }
                else//Hij
                {
                    if(knowledge->get_range(1, sel_hloc)==2)
                    {
                        //printf("unknown ");
                        value=Lower_bound[j]+Interval[j]*RandomFloat();
                    }
                    else if(knowledge->get_range(1, sel_hloc)==1)
                    {
                        //printf("positive ");
                        //PAUSE
                        value=0+Upper_bound[j]*RandomFloat();
                    }
                    else if(knowledge->get_range(1, sel_hloc)==-1)
                    {
                        //printf("negative ");
                        //PAUSE
                        value=Lower_bound[j]+(0-Lower_bound[j])*RandomFloat();
                    }
                    //printf("Hij: %f\n",value);
                }
            }
            else//alpha,,beta,deltat
            {
                value=Lower_bound[j]+Interval[j]*RandomFloat();
            }
            *p1 = value;
            //printf("%f ",value);
        }//end of one chromosome
        //printf("\n");
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
    
    for(i=0;i<NO_POP;i++)
    {
        if(obj->EVAL[i] < obj->EVAL[BestOne])
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
    
    // sort by value of fitness function
    for(i=0;i<NO_POP-1;i++)
    {
        for(j=i+1; j<NO_POP;j++)
        {
            if(obj->EVAL[SortIndex[j]] < obj->EVAL[SortIndex[i]])
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
    for(i=0;i<NumCrossover;i+=2)
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
void IGA::Mutation()
{
    //printf("start mutation\n");
    int i;//nGene;
    int m_pop,m_param,point;
    int g_index;
    int h_index;
    int fix_loc;
    int gij_index;
    int hij_index;
    double gij,hij;
    bool loc_same = false;
    //double *mutant;
    //mutant = new double [PARAM_NUM];
    //memset(mutant, 0, PARAM_NUM);
    double value, oldvalue;
    
    if(!knowledge->get_noregulation())
    {
        //============g_location====================================//
        for(i=0;i<m_size2;i++)
        {
            //printf("Mutate Location");
            do{
                m_pop = RandomInt(NO_POP);
            }while(m_pop == BestOne);
            //printf("m_pop=%d BestOne=%d\n",m_pop,BestOne);
            do{
                m_param = RandomInt(2*NO_g);
            }while(m_param%2 == 1);
            //printf("mutant Gij location=%d\n",m_param);
            //PAUSE
            g_index = m_param;
            oldvalue = *(Chromosomes+PARAM_NUM*m_pop+g_index);
            fix_loc = static_cast<int>(MutateOneValue(oldvalue,g_index,0));
            *(Chromosomes+PARAM_NUM*m_pop+g_index) = fix_loc;
            if(oldvalue == fix_loc)
                loc_same = true;
            if(!knowledge->get_fixg() && !loc_same)
            {
                //printf("loc change, gij have to change!\n");
                //PAUSE
                gij_index = g_index + 1;
                gij = *(Chromosomes + PARAM_NUM*m_pop + gij_index);
                *(Chromosomes+PARAM_NUM*m_pop+gij_index) = MutateOneValue(gij, gij_index, fix_loc);
                loc_same = false;
            }
        }
        //=============Gij========================================//
        for(i=0;i<m_size1;i++)
        {
            do{
                m_pop = RandomInt(NO_POP);
            }while(m_pop == BestOne);
            do{
                m_param = RandomInt(2*NO_g);
            }while(m_param%2 == 0);
            //printf("mutant Gij\n");
            //printf("mutant index=%d\n",m_param);
            //PAUSE
            gij_index = m_param;
            g_index = static_cast<int>(*(Chromosomes+PARAM_NUM*m_pop+(gij_index-1)));
            //printf("beta_loc=%d\n",beta_loc);
            gij = *(Chromosomes+PARAM_NUM*m_pop+gij_index);
            *(Chromosomes+PARAM_NUM*m_pop+gij_index)=MutateOneValue(gij,gij_index,g_index);
        }
        //printf("finish Gij\n");
    }
    
    
    
    //============h_location=================================//
    for(i=0;i<m_size2;i++)
    {
        //printf("mutant location\n");
        do{
            m_pop = RandomInt(NO_POP);
        }while(m_pop == BestOne);
        do{
            m_param = RandomInt(2*NO_h);
        }while(m_param%2 == 1);
        h_index = h_start + m_param;
        oldvalue = *(Chromosomes+PARAM_NUM*m_pop+h_index);
        fix_loc = static_cast<int>(MutateOneValue(oldvalue,h_index,0));
        *(Chromosomes+PARAM_NUM*m_pop+h_index) = fix_loc;
        if(oldvalue == fix_loc)
            loc_same = true;
        if(!knowledge->get_fixh() && !loc_same)
        {
            hij_index = h_index + 1;
            hij = *(Chromosomes + PARAM_NUM*m_pop + hij_index);
            *(Chromosomes+PARAM_NUM*m_pop+hij_index) = MutateOneValue(hij, hij_index, fix_loc);
            loc_same = false;
        }
    }
    //============hij=================================//
    for(i=0;i<m_size1;i++)
    {
        //printf("mutant location\n");
        do{
            m_pop = RandomInt(NO_POP);
        }while(m_pop == BestOne);
        do{
            m_param = RandomInt(2*NO_h);
        }while(m_param%2 == 0);
        hij_index = h_start + m_param;
        h_index = static_cast<int>(*(Chromosomes+PARAM_NUM*m_pop+(hij_index-1)));
        hij = *(Chromosomes+PARAM_NUM*m_pop+hij_index);
        *(Chromosomes+PARAM_NUM*m_pop+hij_index)=MutateOneValue(hij,hij_index,h_index);
    }
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
    int old_loc, tmp_loc,flag,x;
    //printf("index=%d\n",index);
    //PAUSE
    value = 0;
    if(index < h_start)
    {
        if(index % 2 == 0)//g_Location
        {
            if(knowledge->get_fixg())
                value = oldValue;
            else
            {
                old_loc = (int)oldValue;
                if(old_loc == neg_zero){old_loc = 0;}
                else{old_loc = abs(old_loc);}
                if(knowledge->get_choose(0,old_loc)==1)
                {
                    value = old_loc;
                    //printf("fix_reg\n");
                    //PAUSE
                }
                else
                {
                    
                    
                    x = RandomInt(gloc_possible);
                    tmp_loc = sortloc[0][x];
                    /*
                    do{
                        tmp_loc = RandomInt(NumTFs);
                    }while(knowledge->get_LocMask(0,tmp_loc)==1);
                    */
                    flag=rand() % 2;
                    if(flag == 0)//negative location
                    {
                        if(tmp_loc == 0){value = neg_zero;}
                        else{value = -tmp_loc;}
                    }
                    else
                    {
                        value = tmp_loc;
                    }
                    
                }

            }
            //printf("after g Locaion mutate: %f\n",value);
            
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
        if(index % 2 == 0)//h_Location
        {
            if(knowledge->get_fixg())
                value = oldValue;
            else
            {
                old_loc = (int)oldValue;
                if(old_loc == neg_zero){old_loc = 0;}
                else{old_loc = abs(old_loc);}
                if(knowledge->get_choose(1,old_loc)==1)
                {
                    value = old_loc;
                }
                else
                {
                    x = RandomInt(hloc_possible);
                    tmp_loc = sortloc[1][x];
                    
                    /*
                    do{
                        tmp_loc = RandomInt(NumTFs);
                    }while(knowledge->get_LocMask(1,tmp_loc)==1);
                    */
                    flag=rand() % 2;
                    if(flag == 0)//negative location
                    {
                        if(tmp_loc == 0){value = neg_zero;}
                        else{value = -tmp_loc;}
                    }
                    else
                    {
                        value = tmp_loc;
                    }
                    
                }
                
            }
            //printf("after h Locaion mutate: %f\n",value);
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
    bool loc_same = false;
    int Index[2];
    //int DoMutationIndex[PARAM_NUM];
    int g_loc, h_loc, fix_loc;
    if(obj->EVAL[index_one] <= obj->EVAL[index_two])
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
                
                if(i < h_start)
                {
                    if(!knowledge->get_noregulation())
                    {
                        if(i % 2 == 0)//g_Location
                        {
                            
                            fix_loc = static_cast<int>(MutateOneValue(OA->m_Parents[1][i],i,0));
                            if(fix_loc == OA->m_Parents[1][i])
                                loc_same = true;
                            OA->m_Parents[1][i] = fix_loc;
                            OA->m_OATable[OA->m_Height][i]=OA->m_Parents[1][i];
                            DoMutationIndex[i]=1;
                            if(!knowledge->get_fixg() && !loc_same)
                            {
                                //printf("cx mutate loc and gij\n");
                                //PAUSE
                                OA->m_Parents[1][i+1] = MutateOneValue(OA->m_Parents[1][i+1],i+1,fix_loc);
                                OA->m_OATable[OA->m_Height][i+1]=OA->m_Parents[1][i+1];
                                DoMutationIndex[i+1]=1;
                                loc_same = false;
                            }
                            
                        }
                        else//Gij
                        {
                            if(DoMutationIndex[i] == 0)
                            {
                                g_loc = static_cast<int>(OA->m_Parents[1][i-1]);
                                OA->m_Parents[1][i]=MutateOneValue(OA->m_Parents[1][i],i,g_loc);
                                OA->m_OATable[OA->m_Height][i]=OA->m_Parents[1][i];
                                DoMutationIndex[i]=1;
                            }
                        }
                        reCom=true;
                    }
                }
                else if(i>=h_start && i<alpha_index)
                {
                    if(i % 2 == 0)//h_Location
                    {
                        fix_loc = static_cast<int>(MutateOneValue(OA->m_Parents[1][i],i,0));
                        if(fix_loc == OA->m_Parents[1][i])
                            loc_same = true;
                        OA->m_Parents[1][i] = fix_loc;
                        OA->m_OATable[OA->m_Height][i]=OA->m_Parents[1][i];
                        DoMutationIndex[i]=1;
                        if(!knowledge->get_fixh() && !loc_same)
                        {
                            OA->m_Parents[1][i+1] = MutateOneValue(OA->m_Parents[1][i+1],i+1,fix_loc);
                            OA->m_OATable[OA->m_Height][i+1]=OA->m_Parents[1][i+1];
                            DoMutationIndex[i+1]=1;
                            loc_same = false;
                        }

                    }
                    else//Hij
                    {
                        if(DoMutationIndex[i] == 0)
                        {
                            h_loc = static_cast<int>(OA->m_Parents[1][i-1]);
                            OA->m_Parents[1][i]=MutateOneValue(OA->m_Parents[1][i],i,h_loc);
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
