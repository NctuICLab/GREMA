#include "OrthogonalArray.h"
#define			Debug_mode 0
//---------------------------------------------------------------------------------------
COrthogonalArray::COrthogonalArray(Objfunction *objection)
{
	m_OA = NULL;
    m_CutPoint = NULL;
    m_OA_MSE = NULL;
    m_Index = NULL;
    m_flag = NULL;
    GeneValueIndex = NULL;
    m_Parents = NULL;
    m_Width = 0;
    m_Height = 0;
    m_Level = 0;
    m_OAMSEHeight = 0;
    m_OAMSEWidth = 0;
    m_NumUsedParam = 0;
    m_NumParameters = 0;
    m_Chromosomes = 0;
    m_NumSection = 0;
    m_NumSecElements = 0;
    obj = objection;
    //printf("finish init OA\n");
}
//---------------------------------------------------------------------------------------
COrthogonalArray::~COrthogonalArray()
{
	int i;
    if(m_OA){delete [] m_OA;}
    if(m_CutPoint){delete [] m_CutPoint;}
    if(m_OA_MSE){delete [] m_OA_MSE;}
    if(m_Index){delete [] m_Index;}
    if(m_flag){delete [] m_flag;}
    if(GeneValueIndex){delete [] GeneValueIndex;}
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
    
}
//---------------------------------------------------------------------------------------
void COrthogonalArray::Create(unsigned int param, double *chromosomes,int NumElements)
{
	unsigned int i;
    int usefactors;
	m_Level=2;
    int length = 2;
	m_NumParameters=param;//codenlen
    if(m_NumParameters%length == 1)
    {
        usefactors = m_NumParameters/length;
    }
    else
    {
        usefactors = (int)ceil(m_NumParameters/ (float)length);
    }
    //printf("usefactor=%d\n",usefactors);
    //PAUSE
	m_Chromosomes=chromosomes;
	OrthogonalArray(usefactors);
    if(Debug_mode)
    {
        DisplayOATable();
        PAUSE
        printf("finish OA(factors)\n");
    }
	m_OAMSEHeight=m_Height+m_Level+1;//8+2+1=11
	m_OAMSEWidth=m_NumParameters+1;//11+1=12
    //printf("m_OAMSEHeight=%d m_OAMSEWidth=%d\n",m_OAMSEHeight,m_OAMSEWidth);
    m_NumSection=m_NumParameters/NumElements;
    m_NumSecElements=NumElements;
    m_CurUsedHeight=usefactors+1;
    m_CurUsedWidth=usefactors;
    //printf("m_curHeigh=%d,m_curWidth=%d,m_step=%d\n",m_CurUsedHeight,m_CurUsedWidth,m_Step);
    //PAUSE

    m_OATable = new double* [m_OAMSEHeight];
	for(i=0;i<m_OAMSEHeight;i++)
	{
        m_OATable[i]=new double[m_OAMSEWidth];
        memset(m_OATable[i],0,m_OAMSEWidth);
	}
    m_Parents = new double* [2];
    for(i=0;i<2;i++)
        m_Parents[i] = new double[m_NumParameters];
	m_flag=new int[m_Height+3];
    m_Index=new unsigned int[param];
    GeneValueIndex=new int[param];
    m_CutPoint=new unsigned int[m_Width+1];
    m_CutPoint[0]=0;
    for(i=1;i<usefactors;i++)
        m_CutPoint[i]=m_CutPoint[i-1]+length;
    m_CutPoint[i] = m_NumParameters;
    if(Debug_mode)
    {
        for(i=0;i<=usefactors;i++)
            printf("m_CutPoint[%d]=%d ",i,m_CutPoint[i]);
        printf("\n");
        PAUSE
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
	
    m_OA = new unsigned int [m_Height*m_Width];
	dist  = m_Height;//4
	block = m_Height/2;//2
	m = 1;
	for(i=0;i<m_Width;i++) {
		if(i+1  == m ) 
		{
			for(j=0;j<m_Height;j+=dist) {
				for(k=0;k<block;k++) {
                    *(m_OA+(j+k)*m_Width+i) = 0;
                    *(m_OA+(j+block+k)*m_Width+i) = 1;
				}
			}
			m <<= 1;	/* m = m * 2 */
			dist >>= 1;
			block >>= 1;
		}
		else {
			k = m/2;
			for(j=0;j<m_Height;j++) {
				temp=*(m_OA+j*m_Width+i-k) ^ *(m_OA+j*m_Width+k-1);
				*(m_OA+j*m_Width+i)=temp;
			}
		}
	}

}
//----------------------------------------------------------------------------------------------
void COrthogonalArray::DisplayOATable()
{
	unsigned int i,j;
	unsigned int *p;
	p=m_OA;
	for(j=0;j<m_Height;j++)
	{
		for(i=0;i<m_Width;i++)
		{
			printf("%d ",*p++);
		}
		printf("\n");
	}
}
//-------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------

double COrthogonalArray::EvaluateFitness(double *genes)
{
    double  solution[m_NumParameters];
    double  value=0;
	int i;
	for(i=0;i<m_NumParameters;i++)
		solution[i]=genes[i];
	value=obj->fitness(solution);
	return value;
}
//-------------------------------------------------------------------------------------------------
//								Intelligent Crossover
//-------------------------------------------------------------------------------------------------
void COrthogonalArray::IntelliegentCrossover(unsigned int index_one,unsigned int index_two,unsigned int best)
{
    unsigned int		i, j, k , Level;
	unsigned int		OffspringIndex[2];//children
	double				temp;
	double				m_MinMSEValue;
	unsigned int		*pOA,m_MinMSEFactor,m_MinFactorIndex=0,mm;
    double              testvalue;
    best_index = best;
    for(i=0;i<m_Height+3;i++)
        m_flag[i]=0;
    if(Debug_mode)
    {
        printf("crossover start! index1=%d    index2=%d\n",index_one,index_two);
        printf("before cx,%d\n",index_one);
        for(i=0;i<m_NumParameters;i++)
        {
            testvalue = *(m_Chromosomes+m_NumParameters*index_one+i);
            printf("%f ",testvalue);
        //m_Parents[0][i]=m_OATable[OffspringIndex[0]][i];
        //m_Parents[1][i]=m_OATable[OffspringIndex[1]][i];
        }
        printf("\nbefore cx,%d\n",index_two);
        for(i=0;i<m_NumParameters;i++)
        {
        
            testvalue = *(m_Chromosomes+m_NumParameters*index_two+i);
            printf("%f ",testvalue);
        //m_Parents[0][i]=m_OATable[OffspringIndex[0]][i];
        //m_Parents[1][i]=m_OATable[OffspringIndex[1]][i];
        }
    }
	if(Debug_mode)
	{
		printf("\nParents.............\n");
		PrintSolution(m_OATable[0]);
		PrintSolution(m_OATable[m_Height]);
		printf("..........\n");
        PAUSE
	}
    //printf("m_Step=%d\n",m_Step);
    for(i=1;i<m_Height;i++)
    {
		for(j=0;j<m_CurUsedWidth;j++)
        {
			Level=*(m_OA+i*m_Width+j);
            //printf("%d ",Level);
			for(k=m_CutPoint[j];k<m_CutPoint[j+1];k++)
            {
				m_OATable[i][k] = m_Parents[Level][k];
				//printf("m_CutPoint[%d]:m_OATable[%d][%d]=%f ",j,i,k,m_OATable[i][k]);
            }
            //PAUSE
        }
        //printf("\n");
		
		m_OATable[i][m_NumParameters] = EvaluateFitness(m_OATable[i]);
		if(Debug_mode)
        {
            printf("%d:",i);
			PrintSolution(m_OATable[i]);
        }
    }
	//====================================================================================
	//計算每個factor的MSE
	
    m_MinMSEValue= MAXDOUBLE;
	m_MinMSEFactor=0;
	mm=m_Width;
	for(i=0;i<m_CurUsedWidth;i++) 
	{
		for(j=0;j<m_Level;j++)
			m_MES[j]=0;
		pOA=m_OA+i;
		for(j=0;j<m_Height;j++,pOA+=mm)
		{
			m_MES[*pOA]+=m_OATable[j][m_NumParameters];
            //printf("m_MES[%d]=%f ",*pOA,m_MES[*pOA]);
           
		}
        //printf("\n");
        //printf("m_MES[0]=%f m_MES[1]=%f\n",m_MES[0],m_MES[1]);
        //PAUSE
		if(m_MES[0]<m_MES[1])
			Level=0;
		else
			Level=1;
		temp=fabs(m_MES[0]-m_MES[1]);
		if(temp<m_MinMSEValue)
		{
			m_MinMSEValue=temp;
			m_MinMSEFactor=i;
			m_MinFactorIndex=1-Level;
		}
        //printf("factor=%d MED[%d]=%f m_MinMSEValue=%f m_MinMSEFactor=%d m_MinFactorIndex=%d\n",i+1,i+1,temp,m_MinMSEValue,m_MinMSEFactor,m_MinFactorIndex);
		//  產生二個子解
		for(k=m_CutPoint[i];k<m_CutPoint[i+1];k++)
		{
			m_OATable[m_Height+1][k] = m_Parents[Level][k];
			m_OATable[m_Height+2][k] = m_Parents[Level][k];
			//printf("m_OATable[%d][%d]=%f ",m_Height+1,m_GenesIndex[k],m_OATable[m_Height+1][m_GenesIndex[k]]);
            //printf("m_OATable[%d][%d]=%f\n",m_Height+2,m_GenesIndex[k],m_OATable[m_Height+2][m_GenesIndex[k]]);
		}

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
	for(k=m_CutPoint[m_MinMSEFactor];k<m_CutPoint[m_MinMSEFactor+1];k++)
	{
		m_OATable[m_Height+2][k] = m_Parents[m_MinFactorIndex][k];
	}
	//====================================================================================
    for(i=1;i<=2;i++)//算出子解的fitness
	{
		m_OATable[m_Height+i][m_NumParameters] = EvaluateFitness(m_OATable[m_Height+i]);
//		printf("%f ",m_OATable[m_Height+i][m_NumParameters]);
	}
	if(Debug_mode)
	{
		printf("\nOffspring ok........\n");
		PrintSolution(m_OATable[m_Height+1]);
		PrintSolution(m_OATable[m_Height+2]);
        PAUSE
	}
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
			
            if(m_OATable[j][m_NumParameters] <m_OATable[k][m_NumParameters] && m_flag[j]==0)
					k = j;
            //printf("k=%d ",k);
		}
        //printf("\n");
        m_flag[k]=1;
        OffspringIndex[i]=k;
        //printf("OffspringIndex[%d]=%d\n",i,k);
        //PAUSE
	}
	//====================================================================================
	
	for(i=0;i<m_NumParameters;i++)
	{
		*(m_Chromosomes+m_NumParameters*index_one+i) = m_OATable[OffspringIndex[0]][i];
        *(m_Chromosomes+m_NumParameters*index_two+i) = m_OATable[OffspringIndex[1]][i];
        
        //m_Parents[0][i]=m_OATable[OffspringIndex[0]][i];
		//m_Parents[1][i]=m_OATable[OffspringIndex[1]][i];
	}
// ------------------------------------------------------------
	if(index_one==best_index && m_OATable[OffspringIndex[0]][m_NumParameters]>obj->EVAL[best_index])
	{
		printf("  Best One Lost!!error!!");
		//getch();
	}
	obj->EVAL[index_one]=m_OATable[OffspringIndex[0]][m_NumParameters];
	obj->EVAL[index_two]=m_OATable[OffspringIndex[1]][m_NumParameters];
    //printf("index_one:EVAL[%d]=%f index_two:EVAL[%d]=%f\n",index_one,EVAL[index_one],index_two,EVAL[index_two]);
	if(Debug_mode)
	{
		printf("Real Offspring........\n");
		PrintSolution(m_OATable[OffspringIndex[0]]);
		PrintSolution(m_OATable[OffspringIndex[1]]);
	}

	if(obj->EVAL[best_index]>obj->EVAL[index_one])
		best_index=index_one;
    if(obj->EVAL[best_index]>obj->EVAL[index_two])
		best_index=index_two;
    if(Debug_mode)
    {
        printf("after cx,%d\n",index_one);
        for(i=0;i<m_NumParameters;i++)
        {
            testvalue = *(m_Chromosomes+m_NumParameters*index_one+i);
            printf("%f ",testvalue);
        }
        printf("%f\n",obj->EVAL[index_one]);
        printf("after cx,%d\n",index_two);
        for(i=0;i<m_NumParameters;i++)
        {
        
            testvalue = *(m_Chromosomes+m_NumParameters*index_two+i);
            printf("%f ",testvalue);
        }
        printf("%f\n",obj->EVAL[index_two]);
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
		for(j=0;j<m_NumParameters+1;j++)
		{
			printf("%6.4f ",m_OATable[i][j]);
		}
		printf("\n");
	}
	printf("======================MSE================================\n");
	for(i=m_Height;i<m_Height+3;i++)
	{
		for(j=0;j<m_NumParameters+1;j++)
		{
			printf("%6.4f ",m_OATable[i][j]);
		}
		printf("\n");
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
