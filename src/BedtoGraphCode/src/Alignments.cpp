#include "Alignments.h"

//Alignment

    string readName;
    char strandType; //S for same, D for differnt, B for both
    string chr;
    int strand;
    uint32_t start;

string Alignment::getKey()
{
    string key="";
    key=key+readName+"_"+chr+"_"+std::to_string((long long int)strand)+"_"+std::to_string((long long int)start);
    key=key+"_"+cigar+"_"+std::to_string((long long int)editDistance);
    return key;
}

bool Alignment::changeFromString(string str)
{
    vector<string> tmp = my_split(str, ',');
    if(tmp.size()!=6)
        return false;
    setReadName(tmp[0]);
    setStrandType(tmp[1][0]);
    setChr(tmp[2]);
    if(tmp[3][0]=='0')
        setStrand(0);
    else
        setStrand(1);
    string posStr=tmp[3].substr(1,tmp[3].length()-1);
    setStart(atol(posStr.c_str()));
    setCigar(tmp[4]);
    setEditDistance(atoi(tmp[5].c_str()));
    return true;
}

//Alignments

Alignment Alignments::getAlignment(string key)
{
    unordered_map<string,Alignment>::iterator it=alignments.find(key);
    return it->second;
}

void Alignments::insertAlignment(string key, Alignment a)
{
  //alignments.insert(make_pair<string,Alignment>(key,a));
  alignments.insert(make_pair(key,a));
}

bool Alignments::isAlignmentIn(string key)
{
    unordered_map<string,Alignment>::iterator it=alignments.find(key);
    if(it!=alignments.end())
        return true;
    else
        return false;
}

bool Alignments::getAlignmentStrings(string line, vector<string> & strs)
{
    std::vector<std::string> tmp = my_split(line, '\t');    
    if(tmp.size()<4)
        return false;
    strs = my_split(tmp[3], ';');
    if(strs.size()<=0)
        return false;
    return true;
}

