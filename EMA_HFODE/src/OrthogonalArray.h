#ifndef COrthogonalArray_H
#define COrthogonalArray_H
#include "Objfunction.h"
#include "Define.h"
class COrthogonalArray
{
public:
	//COrthogonalArray(Objfunction *objection);
    COrthogonalArray();
	~COrthogonalArray();
	void    Create(unsigned int param);
	void    IntelliegentCrossover(unsigned int index_one,unsigned int index_two,unsigned int best);
    void    SetParent(int, int);
    void    initialOA(Objfunction *objection, double **chrom, double *FitnessValue);
    int     set_bestindex(){return best_index;}
	void	DisplayOATable();
    double  **m_Parents;
    double  **m_OATable;
    unsigned int	m_Width, m_Height;
    
private:
    Objfunction *obj;
    void            OrthogonalArray(int factor_number);
    void            DisplayParameterTableData();
    void            PrintSolution(double *p);
    double          EvaluateFitness(double** genes);
    int             RequiredRows(int f);
    unsigned int    best_index;
    unsigned int    **m_OA;
    //unsigned int	m_CurUsedWidth,m_CurUsedHeight;
    unsigned int    m_Level;
    unsigned int    m_OAMSEHeight, m_OAMSEWidth;
    unsigned int    m_NumParameters;
    unsigned int    *m_CutPoint;
    //float           *m_OA_MSE;
    double          **m_Chromosomes;
    double          m_MES[2];
    double          *PopFitenss;
    //unsigned int	*m_Index;
    int             *m_flag;
    //int             *GeneValueIndex;
    //int             m_NumSection,m_NumSecElements;
    int             block;

};
#endif
