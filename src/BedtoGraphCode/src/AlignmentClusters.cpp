#include "AlignmentClusters.h"

void AlignmentClusters::getClusterKeys(string akey, vector<string> & keys)
{
    auto range = alignmentClusters.equal_range(akey);
    unordered_multimap<string,string>::iterator it;
    for(it=range.first;it!=range.second;it++)
    {
        keys.push_back(it->second);
    }
}

void AlignmentClusters::insert(string akey, string ckey)
{
    vector<string> ckeys;
    getClusterKeys(akey,ckeys);
    bool exist=false;
    for(uint i=0;i<ckeys.size();i++)
    {
        if(ckeys[i].compare(ckey)==0)
            exist=true;
    }
    if(exist==false) 
        alignmentClusters.insert(make_pair(akey,ckey));
}
