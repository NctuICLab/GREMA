#include "OrthogonalArray.h"
#define			Debug_mode 0
//---------------------------------------------------------------------------------------
//COrthogonalArray::COrthogonalArray(Objfunction *objection)
COrthogonalArray::COrthogonalArray()
{
	m_OA = NULL;
    m_NumSlots = 0;
    m_flag = NULL;
    m_Parents = NULL;
    m_OATable = NULL;
    m_Width = 0;
    m_Height = 0;
    m_Level = 0;
    m_OAMSEHeight = 0;
    m_OAMSEWidth = 0;
    m_NumParameters = 0;
    m_Chromosomes = 0;
    //obj = objection;
    best_index =0;
    block = 0;
    obj = NULL;
    PopFitenss = NULL;
    //printf("finish init OA\n");
}
//---------------------------------------------------------------------------------------
COrthogonalArray::~COrthogonalArray()
{
	int i;
    
    if(m_OA)
    {
        for(i=0;i<m_Height;i++)
        {
            delete [] m_OA[i];
        }
        delete [] m_OA;
    }
    
    if(m_flag){
        delete [] m_flag;
    }
    if(m_OATable)
    {
        for(i=0;i<m_OAMSEHeight;i++)
            delete [] m_OATable[i];
        delete [] m_OATable;
    }
    
    if(m_Parents)
    {
        for(i=0;i<2;i++)
            delete [] m_Parents[i];
        delete [] m_Parents;
    }
    /*
    if(obj)
    {
        delete obj;
    }
    
    if(m_Chromosomes)
    {
        delete m_Chromosomes;
    }
    if(PopFitenss)
    {
        delete PopFitenss;
    }
    */
}
//---------------------------------------------------------------------------------------
void COrthogonalArray::initialOA(Objfunction *objection, double **chrom, double *FitnessValue)
{
    obj = objection;
    m_Chromosomes = chrom;
    PopFitenss = FitnessValue;
    
}
//---------------------------------------------------------------------------------------
void COrthogonalArray::Create(unsigned int param, int nslots)
{
	//obj = objection;
    unsigned int i;
    int usefactors;
    //i = 0;
    //usefactors = 0;
	m_Level=2;
	m_NumParameters = param;//NumEncodParam
    m_NumSlots = nslots;
    if(m_NumParameters != (unsigned int)(3*nslots+4))
    {
        fprintf(stderr, "COrthogonalArray::Create: chromosome length %u does not match %d slots\n", param, nslots);
        exit(1);
    }
    //one factor (crossover block) per regulator slot + one for b/TransMax/degrade
    m_Blocks.clear();
    m_BlockSlot.clear();
    for(int s=0;s<nslots;s++)
    {
        std::vector<unsigned int> b;
        b.push_back(2*s);//loc
        b.push_back(2*s+1);//n
        b.push_back(2*nslots+s);//k
        m_Blocks.push_back(b);
        m_BlockSlot.push_back(s);
    }
    {
        std::vector<unsigned int> b;
        b.push_back(3*nslots);//b
        b.push_back(3*nslots+1);//TransMax
        b.push_back(3*nslots+2);//degrade
        m_Blocks.push_back(b);
        m_BlockSlot.push_back(-1);
    }
    usefactors = (int)m_Blocks.size();
	OrthogonalArray(usefactors);
    if(Debug_mode)
    {
        DisplayOATable();
        printf("finish OA(factors)\n");
    }
	m_OAMSEHeight = m_Height+m_Level+1;//8+2+1=11
	m_OAMSEWidth = m_NumParameters+1;//NumParameters + fitness value
    //printf("m_OAMSEHeight=%d m_OAMSEWidth=%d\n",m_OAMSEHeight,m_OAMSEWidth);
    m_OATable = new double* [m_OAMSEHeight];
	for(i=0;i<m_OAMSEHeight;i++)
	{
        m_OATable[i]=new double[m_OAMSEWidth];
        memset(m_OATable[i],0,m_OAMSEWidth*sizeof(double));
	}
    m_Parents = new double* [2];
    for(i=0;i<2;i++)
    {
        m_Parents[i] = new double[m_NumParameters];
        memset(m_Parents[i],0,m_NumParameters*sizeof(double));
    }
	m_flag=new int[m_Height+3];
    if(Debug_mode)
    {
        printf("usefactor=%d\n",usefactors);
        printf("m_Width=%d\n",m_Width);
        printf("m_Height=%d\n",m_Height);
        printf("m_OAMSEHeight=%d m_OAMSEWidth=%d\n",m_OAMSEHeight,m_OAMSEWidth);
        PAUSE;
    }
}
//Copy crossover block `factor` of parent `parent` into `row`.  For a slot block this is
//the slot's loc, n and k genes plus bit `slot` of the mask gene, so the slot stays a unit.
void COrthogonalArray::AssembleBlock(double *row, unsigned int factor, unsigned int parent)
{
    const std::vector<unsigned int> &b = m_Blocks[factor];
    for(size_t x=0;x<b.size();x++)
        row[b[x]] = m_Parents[parent][b[x]];
    int slot = m_BlockSlot[factor];
    if(slot >= 0)
    {
        unsigned long long bit = 1ULL << slot;
        unsigned long long pm = (unsigned long long)m_Parents[parent][m_NumParameters-1];
        unsigned long long cur = (unsigned long long)row[m_NumParameters-1];
        if(pm & bit) cur |= bit; else cur &= ~bit;
        row[m_NumParameters-1] = (double)cur;
    }
}
//An individual with no active slot is not allowed (the GA never creates one): if the
//recombined mask is empty, activate the first slot.
void COrthogonalArray::RepairMask(double *row)
{
    if(m_NumSlots <= 0) return;
    unsigned long long all = (1ULL << m_NumSlots) - 1;
    unsigned long long cur = (unsigned long long)row[m_NumParameters-1];
    if((cur & all) == 0)
        row[m_NumParameters-1] = (double)(cur | 1ULL);
}
void COrthogonalArray::SetParent(int index_one, int index_two)
{
    for(int i=0; i<m_NumParameters; i++){
        m_Parents[0][i] = m_Chromosomes[index_one][i];
        m_Parents[1][i] = m_Chromosomes[index_two][i];
    }
    
}
//---------------------------------------------------------------------------------------
//							The function creates the OA of any levels
//---------------------------------------------------------------------------------------
//			Input parameters
//					level:			the number of levels for each factors
//					factor_number:	the number of factors in this OA
//			Ouput:
//					m_OA (unsigned int *)           The  start pointer of OA
//                  m_Height (usigned int):			restores the number of rows in the OA
//                  m_Width (usgned int):			restores the number of columns in the OA
//			External functions:
//					mem_alloc(the size of the required memory)
//								It would return the start pointer of the allocated memory space.
//                              
//																					Writed by L.S. 
//																					2001.04.28
//----------------------------------------------------------------------------------------------

int COrthogonalArray::RequiredRows(int f)
{
    int row = 1;
    while(f)
    {
        f >>= 1;
        row <<= 1;
    }
    return row;
}

void COrthogonalArray::OrthogonalArray(int factor_number)//OrthogonalArray(7)
{
	int rows;
    rows = RequiredRows(factor_number);
    //printf("rows=%d\n",rows);
    m_Width = factor_number;
    m_Height = rows;
    //printf("m_width=%d, m_height=%d\n",m_Width,m_Height);
    //PAUSE
    int	i,j,m,k;
	int	dist,block;
	int temp;
	
    /*
     m_OATable = new double* [m_OAMSEHeight];
     for(i=0;i<m_OAMSEHeight;i++)
     {
     m_OATable[i]=new double[m_OAMSEWidth];
     memset(m_OATable[i],0,m_OAMSEWidth);
     }

     */
    
    m_OA = new unsigned int* [m_Height];
    for(i=0;i<m_Height;i++)
    {
        m_OA[i] = new unsigned int[m_Width];
        memset(m_OA[i],0,m_Width*sizeof(unsigned int));
    }
    //m_OA = 0;
	dist  = m_Height;//4
	block = m_Height/2;//2
	m = 1;
	for(i=0;i<m_Width;i++)
    {
		if(i+1  == m ) 
		{
			for(j=0;j<m_Height;j+=dist)
            {
				for(k=0;k<block;k++)
                {
                    m_OA[j+k][i] = 0;
                    m_OA[j+block+k][i] = 1;
                    //*(m_OA+(j+k)*m_Width+i) = 0;
                    //*(m_OA+(j+block+k)*m_Width+i) = 1;
				}
			}
			m <<= 1;	/* m = m * 2 */
			dist >>= 1;
			block >>= 1;
		}
		else
        {
			k = m/2;
			for(j=0;j<m_Height;j++)
            {
                temp = m_OA[j][i-k] ^ m_OA[j][k-1];
                m_OA[j][i] = temp;
                //temp=*(m_OA+j*m_Width+i-k) ^ *(m_OA+j*m_Width+k-1);
				//*(m_OA+j*m_Width+i)=temp;
			}
		}
	}

}
//----------------------------------------------------------------------------------------------
void COrthogonalArray::DisplayOATable()
{
	int i,j;
	//unsigned int *p;
	//p=m_OA;
    
	for(j=0;j<m_Height;j++)
	{
		for(i=0;i<m_Width;i++)
		{
			printf("%d ",m_OA[j][i]);
		}
		printf("\n");
	}
}
//-------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------
/*
double COrthogonalArray::EvaluateFitness(double** genes)
{
    double* solution = new double[m_NumParameters];
    double* each_gene = NULL;
    each_gene = *genes;
    double  value=0;
	int i;
	for(i=0;i<m_NumParameters;i++)
		solution[i]=each_gene[i];
	value=obj->fitness(&solution);
	return value;
}
 */
//-------------------------------------------------------------------------------------------------
//								Intelligent Crossover
//-------------------------------------------------------------------------------------------------
void COrthogonalArray::IntelliegentCrossover(unsigned int index_one,unsigned int index_two,unsigned int best)
{
    unsigned int		i, j, k ,x, Level;
	unsigned int		OffspringIndex[2];//children
	double				temp;
	double				m_MinMSEValue;
    unsigned int		m_MinMSEFactor,m_MinFactorIndex=0;
    double              testvalue;
    double              indexone_fitness, indextwo_fitness;
    best_index = best;
    //printf("before cx Best fitnessPopfitness[%d]=%e\n",best_index,PopFitenss[best_index]);
    for(i=0;i<m_Height+3;i++)
        m_flag[i]=0;
    
    if(Debug_mode)
    {
        printf("crossover start! index1=%d    index2=%d\n",index_one,index_two);
        printf("before cx,parent%d\n",index_one);
        for(i=0;i<m_NumParameters;i++)
        {
            testvalue = m_Chromosomes[index_one][i];
            //testvalue = *(m_Chromosomes+m_NumParameters*index_one+i);
            printf("%f ",testvalue);
        //m_Parents[0][i]=m_OATable[OffspringIndex[0]][i];
        //m_Parents[1][i]=m_OATable[OffspringIndex[1]][i];
        }
        printf("\nbefore cx,parent%d\n",index_two);
        for(i=0;i<m_NumParameters;i++)
        {
            testvalue = m_Chromosomes[index_two][i];
            //testvalue = *(m_Chromosomes+m_NumParameters*index_two+i);
            printf("%f ",testvalue);
        //m_Parents[0][i]=m_OATable[OffspringIndex[0]][i];
        //m_Parents[1][i]=m_OATable[OffspringIndex[1]][i];
        }
    }
	if(Debug_mode)
	{
		printf("\nParents.............\n");
		PrintSolution(m_OATable[0]);//index_one
		PrintSolution(m_OATable[m_Height]);//index_two
		printf("..........\n");
        PAUSE;
	}
    //printf("m_Width=%d\n",m_Width);
    for(i=1;i<m_Height;i++)
    {
		for(j=0;j<m_Width;j++)
        {
            Level = m_OA[i][j];
            //the whole slot block (loc,n,k,mask bit) comes from one parent
            AssembleBlock(m_OATable[i], j, Level);
        }
        RepairMask(m_OATable[i]);
        //printf("\n");
        //printf("m_OATable[%d]\n",i);
		//PrintSolution(m_OATable[i]);
        m_OATable[i][m_NumParameters] = obj->OAfitenss(m_OATable, i);
        //PrintSolution(m_OATable[i]);
		if(Debug_mode)
        {
            printf("fill the OAtable\n");
            printf("%d:",i);
			PrintSolution(m_OATable[i]);
            PAUSE;
        }
    }
    //printf("finish caculate OAtable\n");
	//==================================================================================
	//Calcualte all MED of all factor
    m_MinMSEValue= MAXDOUBLE;
	m_MinMSEFactor=0;
	//mm=m_Width;
	for(i=0;i<m_Width;i++)
	{
		for(j=0;j<m_Level;j++)
			m_MES[j]=0;
		//pOA=m_OA+i;
		for(j=0;j<m_Height;j++)
		{
            Level  = m_OA[j][i];
            m_MES[Level] += m_OATable[j][m_NumParameters];
            //printf("m_MES[%d]=%f ",*pOA,m_MES[*pOA]);
           
		}
        //printf("\n");
        //printf("m_MES[0]=%f m_MES[1]=%f\n",m_MES[0],m_MES[1]);
        //PAUSE
		if(m_MES[0] < m_MES[1])
			Level=0;
		else
			Level=1;
		temp=fabs(m_MES[0] - m_MES[1]);
		if(temp<m_MinMSEValue)
		{
			m_MinMSEValue=temp;
			m_MinMSEFactor=i;
			m_MinFactorIndex=1-Level;
		}
        //printf("factor=%d MED[%d]=%f m_MinMSEValue=%f m_MinMSEFactor=%d m_MinFactorIndex=%d\n",i+1,i+1,temp,m_MinMSEValue,m_MinMSEFactor,m_MinFactorIndex);
		//  產生二個子解
        if(Debug_mode)
        {
            printf("start calculate the MED\n");
            printf("m_MinMSEFactor=%d\n",m_MinMSEFactor);
            PAUSE;
        }
		AssembleBlock(m_OATable[m_Height+1], i, Level);
		AssembleBlock(m_OATable[m_Height+2], i, Level);

	}
	//====================================================================================
	//			產生第二個子解
    //printf("m_MinMSEValue=%f m_MinMSEFactor=%d m_MinFactorIndex=%d\n",m_MinMSEValue,m_MinMSEFactor,m_MinFactorIndex);
    if(Debug_mode)
    {
        printf("\nOffspring initial........\n");
        PrintSolution(m_OATable[m_Height+1]);
        PrintSolution(m_OATable[m_Height+2]);
        PAUSE
    }
	AssembleBlock(m_OATable[m_Height+2], m_MinMSEFactor, m_MinFactorIndex);
	RepairMask(m_OATable[m_Height+1]);
	RepairMask(m_OATable[m_Height+2]);
	//====================================================================================
    for(i=1;i<=2;i++)//算出子解的fitness
	{
		//PrintSolution(m_OATable[m_Height+i]);
        m_OATable[m_Height+i][m_NumParameters] = obj->OAfitenss(m_OATable, m_Height+i);
//		printf("%f ",m_OATable[m_Height+i][m_NumParameters]);
	}
	if(Debug_mode)
	{
		printf("\nOffspring ok........\n");
		PrintSolution(m_OATable[m_Height+1]);
		PrintSolution(m_OATable[m_Height+2]);
        DisplayParameterTableData();
        PAUSE
	}
    //PrintSolution(m_OATable[m_Height+1]);
    //PrintSolution(m_OATable[m_Height+2]);
    //DisplayParameterTableData();

	//====================================================================================
    // find the best 2 fitness of chromosome from the all OA's exps and the new offspring 2 fitness
	// using selection sort
	if(m_OATable[0][m_NumParameters]==m_OATable[m_Height+1][m_NumParameters])
		m_flag[m_Height+1]=1;
	if(m_OATable[0][m_NumParameters]==m_OATable[m_Height+2][m_NumParameters])
		m_flag[m_Height+2]=1;
    /*
	if(m_Step>1)
		for(i=1;i<m_Height;i+=m_Step)
			m_flag[i]=1;
	*/
    for(i=0;i<2;i++)
    {
		if(m_flag[0]==1)
			k = 1;
		else
			k =0;
        //printf("initial k=%d ",k);
		for(j=0;j<m_Height+3;j++)
		{
			
            if(SafeFit(m_OATable[j][m_NumParameters]) < SafeFit(m_OATable[k][m_NumParameters]) && m_flag[j]==0)
					k = j;
            //printf("k=%d ",k);
		}
        //printf("\n");
        m_flag[k]=1;
        OffspringIndex[i]=k;
        //printf("OffspringIndex[%d]=%d\n",i,k);
        //PAUSE;
	}
	//====================================================================================
	
	for(i=0;i<m_NumParameters;i++)
	{
		m_Chromosomes[index_one][i] = m_OATable[OffspringIndex[0]][i];
    }
    for(i=0;i<m_NumParameters;i++){
        m_Chromosomes[index_two][i] = m_OATable[OffspringIndex[1]][i];
        //printf("chrom[%d][%d]=OATable[%d][%d]\n",index_two,i,OffspringIndex[1],i);
        //*(m_Chromosomes+m_NumParameters*index_one+i) = m_OATable[OffspringIndex[0]][i];
        //*(m_Chromosomes+m_NumParameters*index_two+i) = m_OATable[OffspringIndex[1]][i];
        
        //m_Parents[0][i]=m_OATable[OffspringIndex[0]][i];
		//m_Parents[1][i]=m_OATable[OffspringIndex[1]][i];
	}
    //printf("index_one is OATalbe index=%d\n",OffspringIndex[0]);
    //printf("index_two is OATalbe index=%d\n",OffspringIndex[1]);
    indexone_fitness = m_OATable[OffspringIndex[0]][m_NumParameters];
    indextwo_fitness = m_OATable[OffspringIndex[1]][m_NumParameters];
    PopFitenss[index_one] = indexone_fitness;
    PopFitenss[index_two] = indextwo_fitness;
// ------------------------------------------------------------
	// (The former "Best One Lost" test compared PopFitenss[best_index] with the value
	//  just stored in it, which can never fail, and the cx check re-evaluated both
	//  offspring with obj->fitness() although the OA row they were copied from was
	//  evaluated on the same data; for a NaN fitness that check was always "true" and
	//  aborted the run.  Fitness values are now always finite, see Objfunction.)
	if(Debug_mode)
	{
		if(SafeFit(indexone_fitness) != SafeFit(obj->fitness(index_one)) ||
		   SafeFit(indextwo_fitness) != SafeFit(obj->fitness(index_two)))
		{
			fprintf(stderr, "Crossover: offspring fitness mismatch\n");
			exit(1);
		}
	}

    //printf("index_two fitness = %e\n",obj->fitness(index_two));
    
    
    //printf("index_two Popfitness[%d]=%e\n",index_two,PopFitenss[index_two]);
    
    //obj->set_EVAL(index_one, indexone_fitness);
    //obj->set_EVAL(index_two, indextwo_fitness);
	//obj->EVAL[index_one]=m_OATable[OffspringIndex[0]][m_NumParameters];
	//obj->EVAL[index_two]=m_OATable[OffspringIndex[1]][m_NumParameters];
    //printf("index_one:EVAL[%d]=%f index_two:EVAL[%d]=%f\n",index_one,EVAL[index_one],index_two,EVAL[index_two]);
	if(Debug_mode)
	{
		printf("Real Offspring........\n");
		PrintSolution(m_OATable[OffspringIndex[0]]);
		PrintSolution(m_OATable[OffspringIndex[1]]);
	}
    if(SafeFit(PopFitenss[best_index]) > SafeFit(PopFitenss[index_one]))
        best_index = index_one;
    if(SafeFit(PopFitenss[best_index]) > SafeFit(PopFitenss[index_two]))
        best_index = index_two;
	//if(obj->EVAL[best_index]>obj->EVAL[index_one])
	//	best_index=index_one;
    //if(obj->EVAL[best_index]>obj->EVAL[index_two])
	//	best_index=index_two;
    if(Debug_mode)
    {
        printf("after cx,%d\n",index_one);
        for(i=0;i<m_NumParameters;i++)
        {
            testvalue = m_Chromosomes[index_one][i];
            //testvalue = *(m_Chromosomes+m_NumParameters*index_one+i);
            printf("%f ",testvalue);
        }
        printf("%f\n", PopFitenss[index_one]);
        printf("after cx,%d\n",index_two);
        for(i=0;i<m_NumParameters;i++)
        {
            testvalue = m_Chromosomes[index_two][i];
            //testvalue = *(m_Chromosomes+m_NumParameters*index_two+i);
            printf("%f ",testvalue);
        }
        printf("%f\n", PopFitenss[index_two]);
    }

    //printf("crossover over!\n");
}
//-------------------------------------------------------------------------------------------------

void COrthogonalArray::DisplayParameterTableData()
{
	unsigned int i,j;
	printf("Dispaly Inside Parameters of OA .......................\n");
	for(i=0;i<m_Height;i++)
	{
		for(j=0;j<m_NumParameters;j++)
		{
			printf("%6.4f ",m_OATable[i][j]);
		}
		printf("%e\n",m_OATable[i][j]);
	}
	printf("======================MSE================================\n");
	for(i=m_Height;i<m_Height+3;i++)
	{
		for(j=0;j<m_NumParameters;j++)
		{
			printf("%6.4f ",m_OATable[i][j]);
		}
		printf("%e\n",m_OATable[i][j]);
	}
}
//-------------------------------------------------------------------------------------------------
void COrthogonalArray::PrintSolution(double *p)
{
	int j;
	for(j=0;j<m_NumParameters+1;j++)
		printf("%f ",p[j]);
	printf("\n");
}
