#include "IGA.h"
#define	Debug_mode  0
IGA::IGA()
{
    OA = NULL;
    hill = NULL;
    knowledge = NULL;
    obj = NULL;
    //same_interval=false;
    target = NULL;
    Chromosomes = NULL;
    PopFitness = NULL;
    LockMaskValue = NULL;
    ChooseValue = NULL;
    RangeValue = NULL;
    SortIndex = NULL;
    Lower_bound = NULL;
    Upper_bound = NULL;
    Interval = NULL;
    IntervalChangeByTemplate = NULL;
    Loc_tmp = NULL;
    eachpop_loc = NULL;
    sortloc = NULL;
    LockMaskValue = NULL;
    ChooseValue = NULL;
    RangeValue = NULL;
    ModelParameters = NULL;
    //BestValueEachGeneration = NULL;
    //NumFitnessFunction = NULL;
    DoMutationIndex = NULL;
    pc = 0;
    ps = 0;
    pm1 =0;
    pm2 = 0;
    NO_GEN =0;
    NO_POP =0;
    NO_REG =0;
    PARAM_NUM =0;
    NumCrossover =0;
    m_size1 = 0;
    m_size2 = 0;
    s_size = 0;
    BestOne = 0;
    NumVars = 0;
    NumRuns = 0;
    NumTrials = 0;
    NumTFs = 0;
    Num_nochange = 0;
    //beta_start = 0;
    TransMax_index = 0;
    k_start = 0;
    b_index = 0;
    degrade_index = 0;
    temperature = 0;
    SA = 0;
    cx1 = 0;
    cx2 = 0;
    Index = 0;
    runnow = 0;
        //printf("finish ini iga\n");
}
IGA::~IGA()
{
    if(target)
    {
        delete [] target;
    }
    if(SortIndex)
    {
        delete [] SortIndex;
    }
    if(Lower_bound)
    {
        delete [] Lower_bound;
    }
    if(Upper_bound)
    {
        delete [] Upper_bound;
    }
    if(Interval)
    {
        delete [] Interval;
    }
    if(IntervalChangeByTemplate)
    {
        delete [] IntervalChangeByTemplate;
    }
    if(Loc_tmp)
    {
        delete [] Loc_tmp;
    }
    if(eachpop_loc)
    {
        delete [] eachpop_loc;
    }
    if(sortloc)
    {
        delete [] sortloc;
    }
    if(DoMutationIndex)
    {
        delete [] DoMutationIndex;
    }
    fclose(ptr_progress);
}
void IGA::set_parameter(int index,int iteration,int run_mode, COrthogonalArray *OAtable,Model *hill_f,Getknowledge *know,Objfunction *objection,double	**OrigChromosomes, double *fitnessvalue, int *Lock, int *Choose, int *Range, double **ModelP)
{
    OA = OAtable;
    hill = hill_f;
    knowledge = know;
    obj = objection;
    Chromosomes = OrigChromosomes;
    PopFitness = fitnessvalue;
    LockMaskValue = Lock;
    ChooseValue = Choose;
    RangeValue = Range;
    ModelParameters = ModelP;
    char str1[100];
    pc = DEFAULT_PC;
    ps = DEFAULT_PS;
    pm1 = DEFAULT_PM1;
    pm2 = DEFAULT_PM2;
    SA = DEFAULT_SA;
    //PARAM_NUM = NumParam;
    //NumGeneParam = NumOrigParam;//1 gene(2:b + deg)
    //loc(i1)|n(i1)|loc(i2)|n(i2)|ki1|ki2|b|TransMax|deg|mask|
    //=============2I============|===I===|1+1+1+1 =>3I+4
    //beta_start = 2*NO_REG;
    k_start = 2*NO_REG;
    b_index = 3*NO_REG;
    TransMax_index = 3*NO_REG+1;
    degrade_index = 3*NO_REG+2;
    //printf("k_start: %d, b_index: %d, TransMax_index: %d, degrade_index: %d",k_start,b_index,TransMax_index,degrade_index);PAUSE;
    NumCrossover = (int)((float)pc*NO_POP);
    m_size1 = (int)(NO_POP * PARAM_NUM * pm1);
    m_size2=(int)(NO_POP * PARAM_NUM * pm2);
    s_size =(int)((float)ps * NO_POP);
    //SA=0.997;//cooling rate
    Index = index;
    //runtime = run;
    test_mode = run_mode;
    loc_possible = knowledge->get_Possible();
    sprintf(str1,"progress%d_step%d.txt",Index,iteration);
    ptr_progress=fopen(str1,"w");
    if(Debug_mode)
    {
        printf("pc=%f\n",pc);
        printf("NO_REG=%d\n",NO_REG);
        printf("PARAM_NUM=%d\n",PARAM_NUM);
        printf("GeneIndex=%d\n",Index);
        printf("NumGeneParam=%d\n",NumGeneParam);
        printf("Generation=%d\n",NO_GEN);
        PAUSE;
        //printf("loc_possible=%d\n",loc_possible);
        //test Chrom
        /*
        for(int i=0; i<NO_POP; i++){
            for(int j=0; j<PARAM_NUM; j++){
                Chromosomes[i][j] = 2;
                //printf("%f ",Chromosomes[i][j]);
            }
            //printf("\n");
        }
        */
    }
}
void IGA::initialIGA()
{
    target = new unsigned int[NumCrossover];
    SortIndex = new int[NO_POP];
    Lower_bound = new double[PARAM_NUM];
    Upper_bound = new double[PARAM_NUM];
    Interval = new double[PARAM_NUM];
    IntervalChangeByTemplate = new double[PARAM_NUM];
    
    Loc_tmp = new int[NumVars];//initial selected location
    eachpop_loc = new int[NumVars];
    sortloc = new int[loc_possible];
    //printf("loc_possible=%d\n",loc_possible);
    //BestValueEachGeneration = new double[NO_GEN];
    //memset(BestValueEachGeneration,0,NO_GEN*sizeof(double));
    //NumFitnessFunction = new unsigned int[NO_GEN];
    //memset(NumFitnessFunction,0,NO_GEN*sizeof(unsigned int));
        
    DoMutationIndex =   new int[PARAM_NUM];
    OA->Create(PARAM_NUM);
    //OA->DisplayOATable();
}
void IGA::RangeVariable()
{
    int i=0;
    ///==============Location & nij============///
    for(i=0;i<k_start;i++)//
    {
        if(i%2==0)//location
        {
            Lower_bound[i]=0;
            Upper_bound[i]=NumTFs-1;//5-1=4
        }
        else//nij
        {
            //Lower_bound[i] = hill->beta_downvalue();
            //Upper_bound[i] = hill->beta_upvalue();
            Lower_bound[i] = hill->n_downvalue();
            Upper_bound[i] = hill->n_upvalue();
            Interval[i]=Upper_bound[i]-Lower_bound[i];
        }
    }
    /*
    //===============betaij===========//
    for(i=beta_start;i<k_start;i++){
        Lower_bound[i] = hill->beta_downvalue();
        Upper_bound[i] = hill->beta_upvalue();
        Interval[i]=Upper_bound[i]-Lower_bound[i];
    }
     */
    //================k====================//
    for(i=k_start;i<b_index;i++){
        Lower_bound[i] = hill->k_downvalue();
        Upper_bound[i] = hill->k_upvalue();
        Interval[i]=Upper_bound[i]-Lower_bound[i];
    }
    //================b====================//
    i=b_index;
    Lower_bound[i] = hill->b_downvalue();
    Upper_bound[i] = hill->b_upvalue();
    Interval[i] = Upper_bound[i] - Lower_bound[i];
    //==============TransMax===============//
    i=TransMax_index;
    Lower_bound[i] = hill->TransMax_downvalue();
    Upper_bound[i] = hill->TransMax_upvalue();
    Interval[i] = Upper_bound[i] - Lower_bound[i];
    //==============degrade================//
    i=degrade_index;
    Lower_bound[i] = hill->degrade_downvalue();
    Upper_bound[i] = hill->degrade_upvalue();
    Interval[i]=Upper_bound[i]-Lower_bound[i];
    //========Mask============///
    i=PARAM_NUM-1;
    Lower_bound[i]=1;
    Upper_bound[i]= (1 << NO_REG)-1;
    /*
    //========deltat============///
    i=PARAM_NUM-1;
    //printf("i=%d\n",i);
    Lower_bound[i] = hill->deltat_downvalue();
    Upper_bound[i] = hill->deltat_upvalue();
    Interval[i]    = Upper_bound[i]-Lower_bound[i];
    if(Interval[i] == 0){same_interval = true;}
     */
    
    //for(int x=0;x<PARAM_NUM;x++){
    //    printf("Low[%d]=%f up[%d]=%f inter[%d]=%f\n",x,Lower_bound[x],x,Upper_bound[x],x,Interval[x]);
    //}
    //PAUSE;
    
}
void IGA::run(int run)
{
    runnow = run;
    int i;
    //int step;
    int NowGeneration;
    int fitness_same;
    double nowfitness=1000000;
    double newfitenss=1000000;
    double difference=100;
    //double *p;
    int print_no = NO_GEN/10;
    //step = obj->get_time_interval();
    //printf("step:%d\n",step);
    //PAUSE
    NowGeneration=0;
    fitness_same=0;
    while((fitness_same <= Num_nochange) && (NowGeneration <= NO_GEN))
    {
        Evaluate();
        newfitenss = PopFitness[BestOne];
        //printf("Generation=%d\tEvaluate\tBest=%d\t%e\n",NowGeneration,BestOne,PopFitness[BestOne]);
        if(NowGeneration>0 && nowfitness < newfitenss)
        {
            fprintf(stderr,"Nowfitenss=%e\tNewfitness=%e\n",nowfitness,newfitenss);
            fprintf(stderr,"Evaluate error...........");
            exit(0);
        }
        Selection();
        //printf("Generation=%d\tSelection\tBest=%d\t%e\n",NowGeneration,BestOne,PopFitness[BestOne]);
        if(NowGeneration>0 && nowfitness < newfitenss)
        {
            fprintf(stderr,"Nowfitenss=%e\tNewfitness=%e\n",nowfitness,newfitenss);
            fprintf(stderr,"Selection error...........");
            exit(0);
        }
        Crossover();
        //printf("Generation=%d\tCrossover\tBest=%d\t%e\n",NowGeneration,BestOne,PopFitness[BestOne]);
        if(NowGeneration>0 && nowfitness < newfitenss)
        {
            fprintf(stderr,"Nowfitenss=%e\tNewfitness=%e\n",nowfitness,newfitenss);
            fprintf(stderr,"Crossover error...........");
            exit(0);
        }
        Mutation();
        //printf("Generation=%d\tmutation\tBest=%d\t%e\n",NowGeneration,BestOne,PopFitness[BestOne]);
        if(NowGeneration>0 && nowfitness < newfitenss)
        {
            fprintf(stderr,"Nowfitenss=%e\tNewfitness=%e\n",nowfitness,newfitenss);
            fprintf(stderr,"Mutation error...........");
            exit(0);
        }
        //Evaluate();
        //printf("best index=%d\tbest fitness=%e\n",BestOne, PopFitness[BestOne]);
        ChangeMutationRange();
        difference = nowfitness - PopFitness[BestOne];
        if(difference==0)
            fitness_same++;
        else
            fitness_same=0;
        nowfitness = PopFitness[BestOne];
        if(NowGeneration % print_no == 0)
        {
            fprintf(ptr_progress,"%4d %4d %4d %6d %e %1.3f\n",runnow,fitness_same,NowGeneration, obj->function_cost, PopFitness[BestOne], obj->best_cc);
            fflush(ptr_progress);
            if(test_mode)
            {
                fprintf(stderr,"%4d %4d %4d %6d %e %1.3f\n",runnow,fitness_same,NowGeneration,obj->function_cost,PopFitness[BestOne],obj->best_cc);
                obj->Decode(BestOne);
                //obj->printCurrentModelParameters(BestOne);
            }
        }

        if(NowGeneration>0 && nowfitness < PopFitness[BestOne])
        {
            fprintf(stderr,"IGA error...........");
            exit(0);
        }
        //printf("finish the Generation%d best index=%d\tbest fitness=%e\n",NowGeneration, BestOne, PopFitness[BestOne]);
        NowGeneration++;
    }
    
    //p=Chromosomes+BestOne*PARAM_NUM;
    obj->Decode(BestOne);//reflash to the best chromosome
    for(i=0;i<NumGeneParam;i++)
        ModelParameters[runnow][i] =obj->GeneParamValues[i];
    ModelParameters[runnow][NumGeneParam] = PopFitness[BestOne];
}
void IGA::initialpop()
{
    //printf("initial start!!\n");
    if(knowledge->get_noregulation())
    {
        for(int i=0; i<NO_POP; i++){
            for(int j=0; j<PARAM_NUM; j++){
                Chromosomes[i][j] = 0;
                //printf("%f ",Chromosomes[i][j]);
            }
            //printf("\n");
        }
        
    }else{
        int i,j,x;
        int Nij;
        //double *p1;
        double value;
        int Nij_possible_zero;
        int Nij_select_loc;
        int Nij_loc_index;
        int tmp_loc;
        int confirm;
        confirm = knowledge->get_Confirm();
        //printf("confirm:%d\n",confirm);PAUSE;
        temperature = 1.0;
        obj->function_cost = 0;
        obj->best_cc=-1;
        value = 0;
        Nij = 0;
        for(int i=0; i<NO_POP; i++){
            for(int j=0; j<PARAM_NUM; j++){
                Chromosomes[i][j] = 0;
                //printf("%f ",Chromosomes[i][j]);
            }
            //printf("\n");
        }
        
        for(x=0;x<NumTFs;x++)
        {
            //printf("LockMaskValue[%d]=%d\n",x,LockMaskValue[x]);
            if(LockMaskValue[x] == 0)
            {
                
                sortloc[Nij]=x;
                //printf("beta=%d\n",beta);
                Nij++;
            }
        }
        //PAUSE;
        if(Debug_mode)
        {
            printf("possible locations:\n");
            for(x=0;x<loc_possible;x++)
            {
                printf("%d:%d ",x,sortloc[x]);
            }
            printf("\n");
            PAUSE;
        }
        
        for(x=0;x<PARAM_NUM;x++){IntervalChangeByTemplate[x]=Interval[x];}
        //p1=Chromosomes;
        for(i=0;i<NO_POP;i++)
        {
            //printf("No%d: ",i);
            Nij_select_loc=0;
            Nij_loc_index=0;
            for(x=0;x<NumTFs;x++){
                Loc_tmp[x] = LockMaskValue[x];}
            for(x=0;x<NumTFs;x++){
                eachpop_loc[x] = LockMaskValue[x];
            }
            
            if(!(knowledge->get_fix()))
            {
                if(Debug_mode)
                {
                    printf("location original:\n");
                    for(x=0;x<NumVars;x++)
                        printf("%d:%d ",x,eachpop_loc[x]);
                }
                Nij_possible_zero = RandomInt(NO_REG - confirm + 1);
                //printf("\n Nij number of possible zero:%d\n",Nij_possible_zero);
                //PAUSE
                for(x=0;x<Nij_possible_zero;x++)
                {
                    //tmp_loc=0;
                    do{
                        tmp_loc =  RandomInt(NumTFs);
                        //printf("tmp_loc=%d\n",tmp_loc);
                    }while(eachpop_loc[tmp_loc]==1||ChooseValue[tmp_loc]==1);
                    eachpop_loc[tmp_loc]=1;
                    //printf("%d:mutant loc:%d\n",k,tmp_loc);
                }
                
                if(Debug_mode)
                {
                    printf("\n After mutate Location:\n");
                    for(x=0;x<NumVars;x++)
                        printf("%d:%d ",x,eachpop_loc[x]);
                    printf("\n");
                    PAUSE;
                }
            }
            
            //printf("%d:",i);
            for(j=0;j<PARAM_NUM;j++)
            {
                
                if(j < k_start) //location and nij
                {
                    if(j%2 == 0)//location
                    {
                        if(knowledge->get_fix())
                        {
                            Nij_select_loc = sortloc[Nij_loc_index];
                            value = sortloc[Nij_loc_index];
                        }
                        else
                        {
                            if(knowledge->get_idTFs())
                            {
                                tmp_loc = sortloc[Nij_loc_index];
                                //printf("\n%d %d\n",gloc_index,tmp_loc);
                                //PAUSE
                                Loc_tmp[tmp_loc]=1;
                                Nij_select_loc = tmp_loc;//select location
                                //value = tmp_loc;
                                
                                if(eachpop_loc[Nij_select_loc] == 1)//mutant to negative
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
                                }while(Loc_tmp[tmp_loc]==1);
                                Loc_tmp[tmp_loc]=1;
                                Nij_select_loc = tmp_loc;
                                value = tmp_loc;
                                
                            }
                            
                        }
                        
                        //printf("Location: %f, Nij_select_loc: %d\n",value, Nij_select_loc);
                        //PAUSE
                        Nij_loc_index++;
                    }
                    else//Nij
                    {
                        if(RangeValue[Nij_select_loc] == 2)
                        {
                            value=Lower_bound[j]+Interval[j]*RandomFloat();
                            //value = Lower_bound[j] + RandomInt(Interval[j]);
                        }
                        else if(RangeValue[Nij_select_loc] == 1)
                        {
                            value=0+Upper_bound[j]*RandomFloat();
                            //value = 0 + RandomInt(Upper_bound[j]);
                        }
                        else if(RangeValue[Nij_select_loc] == -1)
                        {
                            value=Lower_bound[j]+(0-Lower_bound[j])*RandomFloat();
                            //value = Lower_bound[j] + RandomInt(0-Lower_bound[j]);
                        }
                        //printf("Gij: %f\n",value);
                    }
                }
                /*
                else if(j>=beta_start && j<k_start)
                {
                    do
                    {
                        value = RandomInt((int)Upper_bound[j]+1);
                    }while(value<Lower_bound[j]);
                }
                 */
                else if(j==PARAM_NUM-1)
                {
                    if(knowledge->get_fix())
                        value = Upper_bound[j];
                    else
                    {
                        MaxMask = (int)Upper_bound[j]+1;
                        do{
                            tmp_loc = rand()%MaxMask;
                        }while(tmp_loc==0);
                        value = tmp_loc | knowledge->NumFixedLinks;
                    }
                    
                }
                else//b,beta,k,deg
                {
                    value=Lower_bound[j]+Interval[j]*RandomFloat();
                }
                Chromosomes[i][j] = value;
                //printf("%f ",Chromosomes[i][j]);
            }//end of one chromosome
            //printf("\n");
        }//end of POPSIZE
        //printf("finish initial\n");
        //PAUSE
    
    }
    
}
void IGA::ChangeMutationRange()
{
    int i;
    for(i=0;i<PARAM_NUM;i++)
        IntervalChangeByTemplate[i] = Interval[i]* temperature;
    temperature *= SA;
    if(temperature<0.000001)
        temperature=0.01;
}
void IGA::Evaluate()
{
    //printf("Evaluate\n");
    int	i;
    //double fitness = 0;
    for(i=0; i<NO_POP; i++)
    {
        //fitness = obj->fitness(i);
        //printf("fitenss=%f\n",fitness);
        PopFitness[i] = obj->fitness(i);
        //printf("PopFitness[%d]=%f\n",i,PopFitness[i]);
    }
    //PAUSE;
    for(i=0;i<NO_POP;i++)
    {
        if(PopFitness[i] < PopFitness[BestOne])
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
            if(PopFitness[SortIndex[j]] < PopFitness[SortIndex[i]])
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
        {
            Chromosomes[SortIndex[i]][j] = Chromosomes[BestOne][j];
        }
        PopFitness[SortIndex[i]] = PopFitness[BestOne];
        //obj->EVAL[SortIndex[i]]=obj->EVAL[BestOne];
    }
    //printf("Selection best index:%d\n",BestOne);
    //printf("Selection best index:%d\tfitness=%f\n",BestOne,obj->fitness(BestOne));
}

void IGA::Crossover()
{
    //printf("crossover\n");
    int i;
    target[0]=BestOne;//parent1
    for(i=1;i<NumCrossover;i++){
        target[i] = RandomInt(NO_POP);  //parent2
    }
    for(i=0;i<NumCrossover;i+=2)
    {
        while(target[i]==target[i+1])
        {
            target[i+1] = RandomInt(NO_POP);
        }
        //printf("1,target[%d]:%f,target[%d]:%f\n",target[i],obj->EVAL[target[i]],target[i+1],obj->EVAL[target[i+1]]);
        CX_mutation(target[i], target[i+1]);
        //printf("crossover: ==========%d %d=========\n",target[i],target[i+1]);
        //PAUSE
        //printf("2,target[%d]:%f,target[%d]:%f\n",cx1,obj->EVAL[cx1],cx2,obj->EVAL[cx2]);
        //PAUSE
        //printf("CX_mutation best index:%d\tfitness=%e\n",BestOne,obj->fitness(BestOne));
        OA->IntelliegentCrossover(cx1, cx2, BestOne);
        BestOne = OA->set_bestindex();
        //printf("ICX best index:%d\tfitness=%e\n",BestOne,obj->fitness(BestOne));
        
    }
    //printf("crossover best index:%d\tfitness=%f\n",BestOne,obj->fitness(BestOne));
    //printf("crossover best index:%d\n",BestOne);
}
void IGA::Mutation()
{
    //printf("start mutation\n");
    int i;//nGene;
    int m_pop,m_param,point;
    int Nij_Loc_index,Nij_value_index;
    int fix_loc;
    double Nij_value;
    //double *mutant;
    //mutant = new double [PARAM_NUM];
    //memset(mutant, 0, PARAM_NUM);
    double value, oldvalue;
    //============location====================================//
    for(i=0;i<m_size2;i++)
    {
        //printf("Mutate Location");
        do{
            m_pop = RandomInt(NO_POP);
        }while(m_pop == BestOne);
        //printf("m_pop=%d BestOne=%d\n",m_pop,BestOne);
        do{
            m_param = RandomInt(2*NO_REG);
        }while(m_param%2 == 1);
        //printf("mutant Gij location=%d\n",m_param);
        //PAUSE
        Nij_Loc_index = m_param;
        oldvalue = Chromosomes[m_pop][Nij_Loc_index];
        //oldvalue = *(Chromosomes + PARAM_NUM*m_pop + beta_index);
        fix_loc=static_cast<int>(MutateOneValue(oldvalue,Nij_Loc_index,0));
        Chromosomes[m_pop][Nij_Loc_index] = fix_loc;
        //*(Chromosomes + PARAM_NUM*m_pop + beta_index) = fix_loc;
        if(!knowledge->get_fix())
        {
            Nij_value_index = Nij_Loc_index+1;
            Nij_value = Chromosomes[m_pop][Nij_value_index];
            //betaij = *(Chromosomes + PARAM_NUM*m_pop + betaij_index);
            Chromosomes[m_pop][Nij_value_index] = MutateOneValue(Nij_value, Nij_value_index, fix_loc);
            //*(Chromosomes+PARAM_NUM*m_pop+betaij_index) = MutateOneValue(betaij, betaij_index, fix_loc);
        }
        
    }
    //printf("finish location\n");
    //=============Nij========================================//
    for(i=0;i<m_size1;i++)
    {
        do{
            m_pop = RandomInt(NO_POP);
        }while(m_pop == BestOne);
        do{
            m_param = RandomInt(2*NO_REG);
        }while(m_param%2 == 0);
        //printf("mutant Gij\n");
        //printf("mutant index=%d\n",m_param);
        //PAUSE
        Nij_value_index = m_param;
        Nij_Loc_index = static_cast<int>(Chromosomes[m_pop][Nij_value_index-1]);
        //beta_index = static_cast<int>(*(Chromosomes+PARAM_NUM*m_pop+(betaij_index-1)));
        //printf("beta_loc=%d\n",beta_loc);
        Nij_value = Chromosomes[m_pop][Nij_value_index];
        //betaij = *(Chromosomes+PARAM_NUM*m_pop+betaij_index);
        Chromosomes[m_pop][Nij_value_index] = MutateOneValue(Nij_value,Nij_value_index,Nij_Loc_index);
        //*(Chromosomes+PARAM_NUM*m_pop+betaij_index)=MutateOneValue(betaij,betaij_index,beta_index);
    }
    //printf("finish Betaij\n");
    
    //mutate kij,b,TransMax,degrade
    for(i=0;i<m_size1;i++)
    {
        do{
            m_pop = RandomInt(NO_POP);
        }while(m_pop == BestOne);//choose the chromosome to mutation
        point = RandomInt(NO_REG+3);
        m_param = k_start + point;
        value = Chromosomes[m_pop][m_param];
        //value=*(Chromosomes+PARAM_NUM*m_pop+m_param);
        //printf("basel transcription and degrade value=%f\n",value);
        Chromosomes[m_pop][m_param] = MutateOneValue(value,m_param,0);
        //*(Chromosomes+PARAM_NUM*m_pop+m_param)=MutateOneValue(value,m_param,0);
    }
    //mutant mask
    for(i=0;i<m_size2;i++)
    {
        do{
            m_pop = RandomInt(NO_POP);
        }while(m_pop == BestOne);
        m_param = PARAM_NUM - 1;
        value = Chromosomes[m_pop][m_param];
        //value=*(Chromosomes+PARAM_NUM*m_pop+m_param);
        //printf("mask value=%f\n",value);
        Chromosomes[m_pop][m_param] = MutateOneValue(value,m_param,0);
        //*(Chromosomes+PARAM_NUM*m_pop+m_param)=MutateOneValue(value,m_param,0);
    }
    //printf("finish others\n");
    /*
    if(!same_interval)
    {
        for(i=0;i<m_size1;i++)//delta_t
        {
            do{
                m_pop = RandomInt(NO_POP);
            }while(m_pop == BestOne);//choose the chromosome to mutation
            m_param = PARAM_NUM-1;
            value = Chromosomes[m_pop][m_param];
            // value=*(Chromosomes+PARAM_NUM*m_pop+m_param);
            Chromosomes[m_pop][m_param] = MutateOneValue(value,m_param,0);
            // *(Chromosomes+PARAM_NUM*m_pop+m_param)=MutateOneValue(value,m_param,0);
        }
    }
    */
    //printf("mutation best index:%d\tfitness=%f\n",BestOne,obj->fitness(BestOne));
    
    //printf("mutant finish\n");
}
double IGA::MutateOneValue(double oldValue,int index,int fixloc)
{
    
    if(Debug_mode)
    {
        printf("mutate start!\n");
        printf("oldValue=%f\n",oldValue);
        printf("index=%d\n",index);
        printf("fixloc=%d\n",fixloc);
        //PAUSE
    }
    double r, gradient;
    int tmp_loc, old_loc,flag,x;
    double value = 0;
    //printf("index=%d\n",index);
    //PAUSE
    if(index < k_start)
    {
        if(index % 2 == 0)//Location
        {
            if(knowledge->get_fix())
                value = oldValue;
            else
            {
                old_loc = (int)oldValue;
                if(old_loc == neg_zero)
                    old_loc = 0;
                else
                    old_loc = abs(old_loc);
                if(ChooseValue[old_loc] == 1)
                    value = old_loc;
                else
                {
                    x = RandomInt(loc_possible);
                    tmp_loc = sortloc[x];
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
            
        }
        else//Nij
        {
            //if(fixloc == neg_zero){fix = 0;}
            //else{fix = abs(fixloc);}
            //printf("mutate gij\n");
            //printf("before mutant:%f, range[%d]=%d\n",oldValue,fix,range[fix]);
            //PAUSE
            if(fixloc < 0)
                value = oldValue;
            else{
                if(RangeValue[fixloc] == 2)
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
                else if(RangeValue[fixloc] == 1)
                {
                    //printf("positive\n");
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
                else if(RangeValue[fixloc] == -1)
                {
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
                //printf("after mutant beta:%f\n",value);
                //PAUSE
            }
            
        }
        
    }
    /*
    else if(index>=beta_start && index<k_start)
    {
        do
        {
            value = RandomInt((int)Upper_bound[index]+1);
        }while(value < Lower_bound[index]);
    }
     */
    //Mask
    else if(index == PARAM_NUM-1)
    {
        if(knowledge->get_fix())
            value = oldValue;
        else
        {
            do
            {
                x = rand()%MaxMask;
            }while(x == 0);
            value = x | knowledge->NumFixedLinks;
        }
    }
    /*
    else if(index == PARAM_NUM-1)    //deltat
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
    */
    else	//n,k,b,degrade
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
    //printf("CX_mutation\n");
    unsigned int i;
    
    //double value;
    bool	reCom=false;
    int Index[2];
    //int DoMutationIndex[PARAM_NUM];
    int beta_loc, fix_loc;
    if(PopFitness[index_one] <= PopFitness[index_two])
    //if(*obj->get_EVAL(index_one) <= *obj->get_EVAL(index_two))
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
    //printf("cx1=%d, cx2=%d\n",Index[0],Index[1]);
    //printf("%d,%d\n",index_one,index_two);
    //PAUSE
    OA->SetParent(cx1, cx2);
    //OA->m_Parents[0] = Chromosomes[Index[0]];
    //OA->m_Parents[1] = Chromosomes[Index[1]];
    //OA->m_Parents[0]=Chromosomes+Index[0]*PARAM_NUM;
    //OA->m_Parents[1]=Chromosomes+Index[1]*PARAM_NUM;
    if(Debug_mode)
    {
        printf("\nIndex0:%d fitness:%f\n",cx1, PopFitness[cx1]);
        for(i=0;i<PARAM_NUM;i++)
        {
            printf("%f ",OA->m_Parents[0][i]);
        }
        printf("\n");
        printf("Index1:%d fitness:%f\n",cx2, PopFitness[cx2]);
        for(i=0;i<PARAM_NUM;i++)
        {
            printf("%f ",OA->m_Parents[1][i]);
        }
        printf("\n");
        //PAUSE;
    }
    OA->m_OATable[0][PARAM_NUM] = PopFitness[cx1];
    OA->m_OATable[OA->m_Height][PARAM_NUM] = PopFitness[cx2];
    for(i=0;i<PARAM_NUM;i++)
        DoMutationIndex[i]=0;
    //if(same_interval){DoMutationIndex[PARAM_NUM-1]=1;}
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
                
                if(i < k_start)
                {
                    if(i % 2 == 0)//Location
                    {
                        fix_loc = static_cast<int>(MutateOneValue(OA->m_Parents[1][i],i,0));
                        OA->m_Parents[1][i] = fix_loc;
                        OA->m_OATable[OA->m_Height][i]=OA->m_Parents[1][i];
                        DoMutationIndex[i]=1;
                        if(!knowledge->get_fix())
                        {
                            OA->m_Parents[1][i+1] = MutateOneValue(OA->m_Parents[1][i+1],i+1,fix_loc);
                            OA->m_OATable[OA->m_Height][i+1]=OA->m_Parents[1][i+1];
                            DoMutationIndex[i+1]=1;
                        }
                    }
                    else//Nij
                    {
                        if(DoMutationIndex[i]==0)
                        {
                            beta_loc = static_cast<int>(OA->m_Parents[1][i-1]);
                            OA->m_Parents[1][i]=MutateOneValue(OA->m_Parents[1][i],i,beta_loc);
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
        //printf("have reCom\n");
        //PopFitness[cx2] = obj->fitness(cx2);
        //printf("PopFitness[%d]=%e\n",cx2,PopFitness[cx2]);
        //obj->set_EVAL(Index[1], obj->fitness(&OA->m_Parents[1]));
        //OA->m_OATable[OA->m_Height][PARAM_NUM] = PopFitness[cx2];
        //double test_value = 0;
        //test_value = obj->OAfitenss(OA->m_OATable, OA->m_Height);
        OA->m_OATable[OA->m_Height][PARAM_NUM] = obj->OAfitenss(OA->m_OATable, OA->m_Height);
        
    }
    //printf("CX mutation over!!\n");
    //PAUSE
}
