#ifndef COrthogonalArray_H
#define COrthogonalArray_H
#include "Objfunction.h"
#include "Define.h"
#include <vector>
class COrthogonalArray
{
public:
	//COrthogonalArray(Objfunction *objection);
    COrthogonalArray();
	~COrthogonalArray();
	//param: chromosome length; nslots: number of regulator slots (layout
	//loc(s)|n(s) s=0..N-1, k(s) s=0..N-1, b, TransMax, degrade, mask)
	void    Create(unsigned int param, int nslots);
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
    //Crossover blocks.  One OA factor per regulator slot: {loc(s), n(s), k(s)} plus bit s of
    //the mask gene, so a slot always crosses over as a unit and a parent's k/mask can never
    //be attached to a different regulator.  One last factor holds {b, TransMax, degrade}.
    std::vector< std::vector<unsigned int> > m_Blocks;
    std::vector<int>    m_BlockSlot;//slot of each block (-1: not a slot block)
    int             m_NumSlots;
    void            AssembleBlock(double *row, unsigned int factor, unsigned int parent);
    void            RepairMask(double *row);
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
