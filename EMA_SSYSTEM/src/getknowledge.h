#ifndef getknowledge_h
#define getknowledge_h
class Knowledge
{
public:
    Knowledge();
    ~Knowledge();
    int get_gNumconnections(){ return g_NumConnections; }
    int get_hNumconnections(){ return h_NumConnections; }
    int get_LocMask(int i,int j){return LocMask[i][j];}
    int get_range(int i,int j){return range[i][j];}
    int get_choose(int i,int j){return choose[i][j];}
    int get_gConfirm(){return g_regulations;}
    int get_hConfirm(){return h_regulations;}
    int get_gPossible(){return g_possible;}
    int get_hPossible(){return h_possible;}
    bool get_fixg(){return fix_gknowledge;}
    bool get_fixh(){return fix_hknowledge;}
    bool get_idgTF(){return id_gTF;}
    bool get_idhTF(){return id_hTF;}
    bool get_noregulation(){return no_regulations;}
    //Number of gene slots (location + value pairs) of the chromosome and how they are used:
    //  fixed row (no '?')            : NumConnections = possible, every regulator is kept
    //  random subset (idTF == false) : NumConnections = max(max_reg,confirmed) slots; the first
    //                                  'confirmed' slots are dedicated to the '+'/'-' genes, the others
    //                                  hold distinct randomly chosen '?' genes
    //  otherwise (idTF == true)      : one slot per possible gene, the sign of the location code
    //                                  switches the free ones on/off, '+'/'-' genes are always on
    void initial_knowledge();
    void    Readknowledge(int GeneIndex, char* gknow, char* hknow);
    void    RangeVariable();
    int NumVars;
    int NumTFs;
    
private:
    
    void    set_connections(int possible,int confirmed,int max_reg,bool fixed,int &connections,bool &idTF);
    int **range;
    int **choose;
    int **LocMask;
    int	g_NumConnections;
    int g_possible;
    int h_NumConnections;
    int h_possible;
    int g_regulations;
    int h_regulations;//confirm regulations
    int num_noregulate;
    int PARAM_NUM;
    bool   fix_gknowledge;
    bool   fix_hknowledge;
    bool   no_regulations;
    bool   id_gTF;
    bool   id_hTF;
};

#endif
