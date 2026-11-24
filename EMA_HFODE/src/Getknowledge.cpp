#include "Define.h"
#include "Getknowledge.h"
#include <string>
#include <cstdlib>
#define Debug_mode 0
using namespace std;
Getknowledge::Getknowledge()
{
    NumVars = 0;
    NumTFs = 0;
    NumConnections=0;
    regulations=0;//confirm regulations
    possible=0;//including: '?', '+', '-'
    num_noregulate = 0;
    NumFixedLinks = 0;
    fix_knowledge=false;
    no_regulations=false;
    id_TFs=false;
    //range = NULL;
    //choose = NULL;
    //LocMask = NULL;
    //printf("finish ini k\n");
}
Getknowledge::~Getknowledge()
{

}
void Getknowledge::initial_knowledge()
{
    //range = new int [NumVars];
    //memset(range,0,NumVars);
    //choose = new int [NumVars];
    //memset(choose, 0, NumVars);
    //LocMask = new int [NumVars];
    //memset(LocMask, 0, NumVars);
    
}
void Getknowledge::Readknowledge(int GeneIndex, int iter, int *LocMask, int *range, int *choose, char* know)
{
    //printf("Readknowledge\n");
    ifstream fp;
    //int *LocMask = NULL;
    //LocMask = *Loc_info;
    //int *range = NULL;
    //range = *range_info;
    //int *choose = NULL;
    //choose = *choose_info;
    
    string sline;
    char* token;
    char** X;
    //char knowledge[50];
    float val;
    int i, j, set_gene,num;
    iteration = iter;
    set_gene=GeneIndex;
    if(NumVars <= 0 || NumTFs <= 0 || NumTFs > NumVars || set_gene < 0 || set_gene >= NumVars)
    {
        fprintf(stderr, "Readknowledge: bad sizes (num_var=%d, num_tf=%d, gene index=%d)\n",NumVars,NumTFs,set_gene);
        exit(1);
    }
    //printf("%s knowledge file is opende\n",know);
    //sprintf(knowledge, "./data/knowledge_step%d.txt",iteration);
    fp.open(know, ifstream::in);
    if(!fp)
    {
        fprintf(stderr, "%s open failed!!!\n",know);
        exit(1);
    }
    X = new char* [NumVars];
    for(i = 0; i < NumVars; i++)
    {
        X[i] = new char[NumVars];
        memset(X[i], '0', NumVars);//missing entries mean "no connection"
    }
    i=0;
    while(i < NumVars && getline(fp, sline))//no fixed line-length limit
    {
        j=0;
        token = strtok(&sline[0], " \t\r\n");
        while(token != NULL)
        {
            if(j < NumVars)//ignore surplus columns instead of overflowing X[i]
                X[i][j] = *token;
            token = strtok(NULL, " \t\r\n");
            j++;
        }
        if(j > 0)//a blank line is not a row
            i++;
    }
    if(i < NumVars)
        fprintf(stderr, "Warning: %s has only %d rows (num_var=%d), the missing rows are treated as no connection\n",know,i,NumVars);
    num=0;
    for(i=0;i<NumVars;i++)
    {
        for(j=0;j<NumVars;j++)
        {
            //printf("%c ",X[i][j]);
            if((i==j)&&(X[i][j]== '0'))
                num++;
        }
        //printf("\n");
    }
    //printf("num:%d\n",num);
    
    for(j=0 ; j < NumTFs; j++)
    {
        //printf("%c ",X[set_gene][j]);
        char c = X[set_gene][j];
        if(c == '?')
        {
            LocMask[j]=0;
            range[j]=2;//search range include positive, negative, 0
            choose[j]=0;
            possible++;
        }
        else if(c == '+')
        {
            LocMask[j]=0;
            range[j]=1;//search range only positive
            choose[j]=1;
            regulations++;
            possible++;
        }
        else if(c == '-')
        {
            LocMask[j]=0;
            range[j]=-1;//search range only negative
            choose[j]=1;
            regulations++;
            possible++;
        }
        else
        {
            LocMask[j]=1;//no connection
            range[j]=0;
            choose[j]=0;
            num_noregulate++;
        }
        if(Debug_mode)
        {
            printf("LocMask[%d]=%d range[%d]=%d choose[%d]=%d\n",j,LocMask[j],j,range[j],j, choose[j]);
            PAUSE;
        }
    }
    if(regulations > MAX_SLOTS)
    {
        fprintf(stderr, "Too many confirmed regulators (%d), at most %d are supported\n",regulations,MAX_SLOTS);
        exit(1);
    }
    //Confirmed regulators occupy the first `regulations` slots of every chromosome
    //(bits 0..regulations-1 of the slot mask), see IGA::initialpop.
    NumFixedLinks = (1ULL << regulations) - 1;
    for(i = 0; i<NumVars; i++){
        delete [] X[i];
    }
    delete [] X;
    
    /*
    for(i=0;i<NumTFs;i++)
    {
        if(LocMask[i]==0)
            g_possible++;
        else
            num_noregulate++;
    }
     */
    //printf("possible_regulation=%d\n",possible);
    //PAUSE;
    if(num_noregulate==NumTFs)
    {
        no_regulations=true;
        fix_knowledge = true;
    }
    if(possible == regulations)
    {
        fix_knowledge=true;
    }
    //Number of regulator slots (NumConnections) in each chromosome.
    //Confirmed ('+'/'-') regulators ALWAYS get a slot.  When the candidate regulators are
    //unknown (every column is a candidate, or every column but the diagonal) only a random
    //subset of the '?' regulators is searched at a time: the total number of slots is
    //capped at round(possible/3) and this cap only limits the number of ADDITIONAL free
    //slots, i.e. free slots = min(#'?', max(cap - #confirmed, 1)).  So the cap never
    //removes a confirmed regulator, and at least one free slot is kept when a '?' exists.
    if(((NumVars == NumTFs) && (possible == NumVars))//unknown TFs
       || ((num == NumVars)&&(NumVars == NumTFs)&&(possible == (NumVars-1))))//dream challenge,unknown TFs
    {
        int freeall = possible - regulations;
        int freeslots = 0;
        //experiment: GREMA_SLOT_DIV (default 3) divides the number of candidates to get the slot cap
        double slot_div = 3.0;
        if(getenv("GREMA_SLOT_DIV") && atof(getenv("GREMA_SLOT_DIV")) > 0)
            slot_div = atof(getenv("GREMA_SLOT_DIV"));
        val = possible/slot_div;
        int cap = (int)(val + 0.5);
        if(freeall > 0)
        {
            freeslots = cap - regulations;
            if(freeslots < 1) freeslots = 1;
            if(freeslots > freeall) freeslots = freeall;
        }
        NumConnections = regulations + freeslots;
        if(NumConnections > MAX_SLOTS)
        {
            freeslots = MAX_SLOTS - regulations;
            if(freeslots < 1)
            {
                fprintf(stderr, "Too many regulator slots (%d confirmed), at most %d are supported\n",regulations,MAX_SLOTS);
                exit(1);
            }
            NumConnections = MAX_SLOTS;
            fprintf(stderr, "Warning: only %d regulator slots are supported, searching %d free regulators at a time\n",MAX_SLOTS,freeslots);
        }
        id_TFs = false;
    }
    else
    {
        //every '?' regulator has its own slot (no random subset)
        NumConnections = possible;
        if(NumConnections > MAX_SLOTS)
        {
            fprintf(stderr, "Too many candidate regulators (%d), at most %d are supported\n",NumConnections,MAX_SLOTS);
            exit(1);
        }
        id_TFs = true;
    }
    fp.close();
}

