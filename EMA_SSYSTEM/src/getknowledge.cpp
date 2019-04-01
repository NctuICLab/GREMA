#include "stdafx.h"
#include "getknowledge.h"
#define Debug_mode 0
using namespace std;
Knowledge::Knowledge()
{
    NumVars = 0;
    NumTFs = 0;
    g_NumConnections=0;
    h_NumConnections=0;
    g_regulations=0;//confirm regulations
    h_regulations=0;
    g_possible=0;
    h_possible=0;
    num_noregulate=0;
    fix_gknowledge=false;
    fix_hknowledge=false;
    no_regulations=false;
    id_gTF=false;
    id_hTF=false;
    range = NULL;
    choose = NULL;
    LocMask = NULL;
    //printf("finish ini k\n");
}
Knowledge::~Knowledge()
{
    if(LocMask){delete [] LocMask;}
    if(range){delete [] range;}
    if(choose){delete [] choose;}
}
void Knowledge::initial_knowledge()
{
    int i,j;
    range = new int* [2];
    choose = new int* [2];
    LocMask = new int* [2];
    for(i=0;i<2;i++)
    {
        range[i]=new int[NumVars];
        choose[i]=new int[NumVars];
        LocMask[i]=new int[NumVars];
        for(j=0;j<NumVars;j++)
        {
            range[i][j]=0;
            choose[i][j]=0;
            LocMask[i][j]=0;
        }
    }

}
void Knowledge::Readknowledge(int GeneIndex, char* gknow, char* hknow)
{
    //printf("Readknowledge\n");
    fstream fp,fp2;
    char line[BUF_SIZE], tmp[NumVars];
    char* token;
    char** X;
    char** Y;
    //char g_knowledge[50];
    //char h_knowledge[50];
    int i, j, set_gene, max_reg;
    float val;
	val = NumVars/3.0;
	max_reg = (int)(val + 0.5);
    set_gene=GeneIndex;
    //sprintf(g_knowledge,"./data/g_knowledge_step%d.txt",iteration);
    fp.open(gknow,ios::in);
    //sprintf(h_knowledge,"./data/h_knowledge_step%d.txt",iteration);
    fp2.open(hknow, ios::in);
    X = new char* [NumVars];
    for(i = 0; i < NumVars; i++)
        X[i] = new char[NumVars];
    Y = new char* [NumVars];
    for(i = 0; i < NumVars; i++)
        Y[i] = new char[NumVars];
    if(!fp)
    {
        fprintf(stderr, "%s open failed!!!\n",gknow);
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
    if(!fp2)
    {
        fprintf(stderr, "%s open failed!!!\n",hknow);
        exit(0);
    }
    i=0;
    while(fp2.getline(line,sizeof(line),'\n'))
    {
        j=0;
        token = strtok(line, " \t\r\n");
        while(token != NULL)
        {
            Y[i][j] = *token;
            token = strtok(NULL, " \t\r\n");
            j++;
        }
        i++;
    }
    /*
     for(i=0;i<NumVars;i++)
     {
         for(j=0;j<NumVars;j++)
         {
             printf("%c ",X[i][j]);
         }
         printf("\n");
     }
     PAUSE
    */
    
    for(j=0 ; j < NumTFs; j++)
    {
        //printf("%c ",X[set_gene][j]);
        tmp[j] = X[set_gene][j];
        if(tmp[j] == '?')
        {
            LocMask[0][j]=0;
            range[0][j]=2;//search range include positive, negative, 0
            choose[0][j]=0;
            g_possible++;
        }
        else if(tmp[j] == '+')
        {
            LocMask[0][j]=0;
            range[0][j]=1;//search range only positive
            choose[0][j]=1;
            g_regulations++;
            g_possible++;
        }
        else if(tmp[j] == '-')
        {
            LocMask[0][j]=0;
            range[0][j]=-1;//search range only negative
            choose[0][j]=1;
            g_regulations++;
            g_possible++;
        }
        else
        {
            LocMask[0][j]=1;//no connection
            range[0][j]=0;
            choose[0][j]=0;
            num_noregulate++;
        }
        if(Debug_mode)
        {
            printf("LocMask[0][%d]=%d range[0][%d]=%d choose[0][%d]=%d\n",j,LocMask[0][j],j,range[0][j],j, choose[0][j]);
        }
    }
    
    //printf("possible_regulation:%d confirm:%d no_regulations:%d\n",g_NumConnections,g_regulations,num_noregulate);
    //PAUSE
    if(num_noregulate==NumTFs)
    {
        no_regulations=true;
        fix_gknowledge=true;
        //printf("no regulation\n");
    }
    if(g_possible == g_regulations){fix_gknowledge=true;}
    if((NumVars == NumTFs) && (g_possible >= max_reg))
    {
        //val = g_possible/3.0;
        //g_NumConnections = (int)(val + 0.5);
	g_NumConnections = max_reg;
        id_gTF = false;
    //}else if((NumVars == NumTFs) && (g_possible > NumVars/2)){
	//val = NumVars/3.0;
	//g_NumConnections = (int)(val + 0.5);
    }
    else
    {
        g_NumConnections = g_possible;
        id_gTF = true;
    }
    printf("Number of g_NumConnections:%d\n",g_NumConnections);
    //printf("g_confirm=%d\n",g_regulations);
    //PAUSE
    //===============read h_knowledge.txt================================
    for(j=0 ; j < NumVars; j++)
    {
        //printf("%c ",X[set_gene][j]);
        tmp[j] = Y[set_gene][j];
        if(tmp[j] == '?')
        {
            LocMask[1][j]=0;
            range[1][j]=2;//search range include positive, negative, 0
            choose[1][j]=0;
            h_possible++;
        }
        else if(tmp[j] == '+')
        {
            LocMask[1][j]=0;
            range[1][j]=1;//search range only positive
            choose[1][j]=1;
            h_possible++;
            h_regulations++;
        }
        else if(tmp[j] == '-')
        {
            LocMask[1][j]=0;
            range[1][j]=-1;//search range only negative
            choose[1][j]=1;
            h_possible++;
            h_regulations++;
        }
        else
        {
            LocMask[1][j]=1;//no connection
            range[1][j]=0;
            choose[1][j]=0;
        }
        if(Debug_mode)
        {
            printf("LocMask[1][%d]=%d range[1][%d]=%d choose[1][%d]=%d\n",j,LocMask[1][j],j,range[1][j],j, choose[1][j]);
        }
    }
    if(h_possible == h_regulations){fix_hknowledge=true;}
    //PAUSE

    
    
    
    
    if((NumVars == NumTFs) && (h_possible >= max_reg))
    {
        //val = h_possible/3.0;
        //h_NumConnections = (int)(val + 0.5);
	h_NumConnections = max_reg;
        id_hTF = false;
    //}else if((NumVars == NumTFs) && (h_possible > NumVars/2)){
	//	val = NumVars/3.0;
	//	h_NumConnections = (int)(val + 0.5);
	}
    else
    {
        h_NumConnections = h_possible;
        id_hTF = true;
    }
    printf("Number of h_NumConnections:%d\n",h_NumConnections);
    //printf("h_confirm=%d\n",h_regulations);
    fp.close();
    fp2.close();
    //PAUSE
}

