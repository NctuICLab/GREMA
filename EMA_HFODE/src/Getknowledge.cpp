#include "Define.h"
#include "Getknowledge.h"
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
    
    char line[BUF_SIZE], tmp[NumVars];
    char* token;
    char** X;
    //char knowledge[50];
    float val;
    int i, j, set_gene,num;
    int mask=0;
    iteration = iter;
    set_gene=GeneIndex;
    //printf("%s knowledge file is opende\n",know);
    //sprintf(knowledge, "./data/knowledge_step%d.txt",iteration);
    fp.open(know, ifstream::in);
    X = new char* [NumVars];
    for(i = 0; i < NumVars; i++)
        X[i] = new char[NumVars];
    if(!fp)
    {
        fprintf(stderr, "%s open failed!!!\n",know);
        exit(0);
    }
    i=0;
    while(fp.getline(line,sizeof(line),'\n'))
    {
        j=0;
        token = strtok(line, " \t\r\n");
        while(token != NULL)
        {
            X[i][j] = *token;
            token = strtok(NULL, " \t\r\n");
            j++;
        }
        i++;
    }
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
        tmp[j] = X[set_gene][j];
        if(tmp[j] == '?')
        {
            LocMask[j]=0;
            range[j]=2;//search range include positive, negative, 0
            choose[j]=0;
            possible++;
        }
        else if(tmp[j] == '+')
        {
            LocMask[j]=0;
            range[j]=1;//search range only positive
            choose[j]=1;
            regulations++;
            possible++;
            mask=(mask << 1 )+ 1;
        }
        else if(tmp[j] == '-')
        {
            LocMask[j]=0;
            range[j]=-1;//search range only negative
            choose[j]=1;
            regulations++;
            possible++;
            mask=(mask << 1 )+ 1;
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
    NumFixedLinks = mask;
    if(X){
        for(i = 0; i<NumVars; i++){
            delete [] X[i];
        }
        delete [] X;
    }
    
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
    if((NumVars == NumTFs) && (possible == NumVars))//unknown TFs
    {
        val = possible/3.0;
        NumConnections = (int)(val + 0.5);
        id_TFs = false;
    }
    else if((num == NumVars)&&(NumVars == NumTFs)&&(possible == (NumVars-1)))
    {
        //dream challenge,unknown TFs
        val = possible/3.0;
        NumConnections = (int)(val + 0.5);
        id_TFs = false;
    }
    else
    {
	//if(NumVars <= 10){
            NumConnections = possible;
	//}else{
        //   val = possible/2.0;
        //   NumConnections = (int)(val + 0.5);
        //}
        id_TFs = true;
    }
    fp.close();
}

