#ifndef Getknowledge_h
#define Getknowledge_h
class Getknowledge
{
public:
    Getknowledge();
    ~Getknowledge();
    int get_Numconnections(){ return NumConnections; }
    int get_Confirm(){return regulations;}
    int get_Possible(){return possible;}
    bool get_fix(){return fix_knowledge;}
    bool get_noregulation(){return no_regulations;}
    bool get_idTFs(){return id_TFs;}
    void initial_knowledge();
    void    Readknowledge(int GeneIndex, int iter, int *Loc_info, int *range_info, int *choose_info, char* know);
    void    RangeVariable();
    int NumVars;
    int NumTFs;
    unsigned int	NumFixedLinks;
private:
    int	NumConnections;
    int regulations;
    int possible;
    int num_noregulate;
    int PARAM_NUM;
    int iteration;
    bool   fix_knowledge;
    bool   no_regulations;
    bool   id_TFs;
};
#endif
