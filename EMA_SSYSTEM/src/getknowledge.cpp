#include "stdafx.h"
#include "getknowledge.h"
#include <string>
#include <vector>
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
    int i;
    for(i=0;i<2;i++)
    {
        if(LocMask){delete [] LocMask[i];}
        if(range){delete [] range[i];}
        if(choose){delete [] choose[i];}
    }
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
//Read one knowledge matrix (rows = target genes, columns = regulators, one char per entry).
//The file may contain blank lines; every other line must have exactly `cols` entries and there must
//be exactly `rows` of them. Anything else is reported and stops the program (it used to write
//outside of the matrix).
static void read_knowledge_matrix(const char *name, char **M, int rows, int cols)
{
    std::ifstream fp(name);
    std::string line;
    int i = 0, j;
    char *token;
    if(!fp)
    {
        fprintf(stderr, "%s open failed!!!\n", name);
        exit(1);
    }
    while(std::getline(fp, line))
    {
        std::vector<char> buf(line.begin(), line.end());
        buf.push_back('\0');
        j = 0;
        token = strtok(&buf[0], " \t\r\n");
        if(token == NULL)
            continue;//blank line
        if(i >= rows)
        {
            fprintf(stderr, "%s: more than %d rows\n", name, rows);
            exit(1);
        }
        while(token != NULL)
        {
            if(j >= cols)
            {
                fprintf(stderr, "%s: row %d has more than %d entries\n", name, i+1, cols);
                exit(1);
            }
            M[i][j] = *token;
            token = strtok(NULL, " \t\r\n");
            j++;
        }
        if(j != cols)
        {
            fprintf(stderr, "%s: row %d has %d entries, expected %d\n", name, i+1, j, cols);
            exit(1);
        }
        i++;
    }
    if(i != rows)
    {
        fprintf(stderr, "%s: found %d rows, expected %d\n", name, i, rows);
        exit(1);
    }
}
//Decide how many gene slots the chromosome has, see getknowledge.h
void Knowledge::set_connections(int possible,int confirmed,int max_reg,bool fixed,int &connections,bool &idTF)
{
    int limit = (max_reg > confirmed) ? max_reg : confirmed;//max_reg only limits the additional free regulators
    if(fixed)
    {
        connections = possible;//keep every regulator of a fully determined row
        idTF = true;
    }
    else if((NumVars == NumTFs) && (possible >= limit))
    {
        connections = limit;
        idTF = false;
    }
    else
    {
        connections = possible;
        idTF = true;
    }
}
void Knowledge::Readknowledge(int GeneIndex, char* gknow, char* hknow)
{
    char** X;
    char** Y;
    int i, j, max_reg;
    float val;
    val = NumVars/3.0;
    max_reg = (int)(val + 0.5);
    if(GeneIndex < 0 || GeneIndex >= NumVars)
    {
        fprintf(stderr, "gene index %d is out of range (0..%d)\n", GeneIndex, NumVars-1);
        exit(1);
    }
    if(NumTFs < 0 || NumTFs > NumVars)
    {
        fprintf(stderr, "num_tf=%d must be between 0 and num_var=%d\n", NumTFs, NumVars);
        exit(1);
    }
    X = new char* [NumVars];
    Y = new char* [NumVars];
    for(i = 0; i < NumVars; i++)
    {
        X[i] = new char[NumTFs+1];
        Y[i] = new char[NumVars+1];
        memset(X[i], '0', NumTFs+1);
        memset(Y[i], '0', NumVars+1);
    }
    read_knowledge_matrix(gknow, X, NumVars, NumTFs);
    read_knowledge_matrix(hknow, Y, NumVars, NumVars);
    //regulators that are not transcription factors can never regulate
    for(j = NumTFs; j < NumVars; j++)
    {
        LocMask[0][j]=1;
        range[0][j]=0;
        choose[0][j]=0;
    }
    for(j=0 ; j < NumTFs; j++)
    {
        //printf("%c ",X[set_gene][j]);
        char c = X[GeneIndex][j];
        if(c == '?')
        {
            LocMask[0][j]=0;
            range[0][j]=2;//search range include positive, negative, 0
            choose[0][j]=0;
            g_possible++;
        }
        else if(c == '+')
        {
            LocMask[0][j]=0;
            range[0][j]=1;//search range only positive
            choose[0][j]=1;
            g_regulations++;
            g_possible++;
        }
        else if(c == '-')
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
    set_connections(g_possible,g_regulations,max_reg,fix_gknowledge,g_NumConnections,id_gTF);
    if(g_regulations > max_reg)
        printf("g: %d confirmed regulators exceed max_reg=%d, all of them are kept\n",g_regulations,max_reg);
    printf("Number of g_NumConnections:%d\n",g_NumConnections);
    //printf("g_confirm=%d\n",g_regulations);
    //PAUSE
    //===============read h_knowledge.txt================================
    for(j=0 ; j < NumVars; j++)
    {
        //printf("%c ",X[set_gene][j]);
        char c = Y[GeneIndex][j];
        if(c == '?')
        {
            LocMask[1][j]=0;
            range[1][j]=2;//search range include positive, negative, 0
            choose[1][j]=0;
            h_possible++;
        }
        else if(c == '+')
        {
            LocMask[1][j]=0;
            range[1][j]=1;//search range only positive
            choose[1][j]=1;
            h_possible++;
            h_regulations++;
        }
        else if(c == '-')
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

    
    
    
    
    set_connections(h_possible,h_regulations,max_reg,fix_hknowledge,h_NumConnections,id_hTF);
    printf("Number of h_NumConnections:%d\n",h_NumConnections);
    //printf("h_confirm=%d\n",h_regulations);
    for(i = 0; i < NumVars; i++)
    {
        delete [] X[i];
        delete [] Y[i];
    }
    delete [] X;
    delete [] Y;
}

