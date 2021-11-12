#include "Groups.h"

int Member::getMemberType(ClusterGraph & cf,Exons & es, string mkey, int numAlignM, vector<string> mAnnotStr, float p1, float p2, float p3)
{
    vector<string> cAnnotStr=getAnnotStr(es);
    ckey=getCkey();
  
    //numAlignments mAlignM
    int share=cf.getReadsNum(mkey, ckey);
    //int res=0;
    if(numAlignM*p1<numAlignments && share<numAlignments*p2) //large enough and not much shared
        return 2;
        //res=2;
    /*if(res==2)
    {
        int numshare=0;
        vector<string> cAnnotStr=getAnnotStr(es);
        for(uint i=0;i<mAnnotStr.size();i++)
            for(uint j=0;j<cAnnotStr.size();j++)
            {
                if(mAnnotStr[i].compare(cAnnotStr[j])==0)
                    numshare++;
            }
        if(numshare>0)
            return 1;
        else 
            return 2;
    }
    */
    if(share>numAlignM*(1.0-p3)) //close number
        return 1;
    return 0;
}

vector<string> Member::getAnnotStr(Exons & es)
{
    vector<string> res;

    unordered_set<string> tmp;
    for(uint i=0;i<eIds.size();i++)
    {
        Exon e=es.getExon(eIds[i]);
        string tmp1=e.getGName();
        tmp.insert(tmp1);
    }

    unordered_set<string>::iterator it;
    for(it=tmp.begin();it!=tmp.end();it++)
    {
        res.push_back(*it);
    }
    return res; 
}

string Member::getMemberStr(Exons & es,ClusterGraph & cf, string mkey)
{
    string res;
    res+=ckey;
    res+="\t"+to_string(static_cast<long long>(numAlignments));
    res+="\t"+to_string(static_cast<long long>(degree));  
    
    res+="\t"+to_string(static_cast<long long>(cf.getReadsNum(mkey, ckey)));
 
    unordered_set<string> tmp;
    for(uint i=0;i<eIds.size();i++)
    {
        Exon e=es.getExon(eIds[i]);
        string tmp1=e.getGName()+"("+e.getGBioType()+")";
        tmp.insert(tmp1);
    }
    
    unordered_set<string>::iterator it;
    for(it=tmp.begin();it!=tmp.end();it++)
    {
        if(it==tmp.begin())
            res+="\t";
        if(it!=tmp.begin())
            res+=";";
        res+=*it;
    }
    res+="\n"; 
    return res;
}

string Member::getMemberStr2(Exons & es,ClusterGraph & cf, string mkey)
{
    string res;
    res+=ckey;
    res+=";degree:"+to_string(static_cast<long long>(degree));
    res+=";alignments:"+to_string(static_cast<long long>(numAlignments));
    res+=";uniq:"+to_string(static_cast<long long>(numUniqAlignments));;
    if(mkey.compare(ckey)==0)
        res+=";overlap:"+to_string(static_cast<long long>(numAlignments));
    else
        res+=";overlap:"+to_string(static_cast<long long>(cf.getReadsNum(mkey, ckey)));

    unordered_set<string> tmp;
    for(uint i=0;i<eIds.size();i++)
    {
        Exon e=es.getExon(eIds[i]);
        string tmp1=e.getGName()+"("+e.getGBioType()+")";
        tmp.insert(tmp1);
    }

    unordered_set<string>::iterator it;
    for(it=tmp.begin();it!=tmp.end();it++)
    {
        if(it==tmp.begin())
            res+=";";
        if(it!=tmp.begin())
            res+=",";
        res+=*it;
    }
    //res+="\n";
    return res;
}


void Group::setMembers(vector<string> & ckeys, ClusterGraph & cg, Alignments & al, ReadAlignments & ra ,Clusters & c, Exons & es)
{
    for(uint i=0;i<ckeys.size();i++)
    {
        Member mb;
        mb.setCkey(ckeys[i]);
        if(ckeys.size()==1)
            mb.setDegree(0);
        else
            mb.setDegree(cg);
        mb.setNumAlignments(c);
        mb.setNumUniqAlignments(c,al, ra);    

        Cluster a=c.getCluster(ckeys[i]);
        vector<uint> ekeys;
        es.coordinateToExons(a.getChr(),a.getStart(),a.getStart()+a.getLength()-1,a.getStrand(),ekeys);
        mb.setEIds(ekeys);
        mb.setStatus(0);
        members.push_back(mb);
    }
}

bool my_sort_member_func(Member i, Member j)
{
    if(i.getNumAlignments()>j.getNumAlignments())
    {
        return true;
    }
    else if(i.getNumAlignments()==j.getNumAlignments())
    {
        if(i.getDegree()>j.getDegree())
        {
            return true;
        }
        else if(i.getDegree()==j.getDegree())
        {
             if((i.getEIds()).size()>0 && (j.getEIds()).size()==0)
             {
                 return true;
             }
             else
             {
                return false;
             }
        }
        else
        {
            return false;
        }
    }
    else
    {
        return false;
    }
    
}

void Group::sortMembers()
{
    sort(members.begin(),members.end(),my_sort_member_func);
}

int Group::getAllAlignmentSize()
{
    int res=0;
    for(uint i=0;i<members.size();i++)
        res+=members[i].getNumAlignments();
    return res; 
}

string Group::getAllMemStr(Exons & es,ClusterGraph & cf)
{
    string res="";
    res=res+to_string(static_cast<long long>(members.size()))+"\n";
    string mkey = members[0].getCkey();
    for(uint i=0;i<members.size();i++)
        res=res+members[i].getMemberStr(es,cf,mkey);
    return res;
}

string Group::getAllMemAbsorbStr(Exons & es,ClusterGraph & cf)
{
    string res="";
    res=res+to_string(static_cast<long long>(members.size()))+":";
    string mkey = members[0].getCkey();
    for(uint i=0;i<members.size();i++)
    {
        if(i!=0)
            res=res+"|";
        res=res+members[i].getMemberStr2(es,cf,mkey);
    }
    return res;
}

string string_append(string str, string append, string sep)
{
    if(str.compare("")==0 || append.compare("")==0)
        return append;
    else
        return str+sep+append;
}

string get_bed_12_str(string key, Clusters & ct)
{
//cout<<"in get_bed_12_str"<<endl;
    uint32_t start=ct.getCluster(key).getStart()-1;
    uint32_t end=start+ct.getCluster(key).getLength();    
    string startStr=to_string(static_cast<long long>(start));
    string endStr=to_string(static_cast<long long>(end));
//cout<<"in get 1"<<endl;
    string res="";
//cout<<"key="<<key<<endl;    
    res+=ct.getCluster(key).getChr();
//cout<<"in get 1.1"<<endl;
    res+="\t"+startStr;
    res+="\t"+endStr;
    res+="\t"+key;
    res+="\t0";
//cout<<"in get 2"<<endl;
    int intStrand=ct.getCluster(key).getStrand();
    string strStrand=".";
    if(intStrand==0)
        strStrand="+";
    if(intStrand==1)
        strStrand="-";
    res+="\t"+strStrand;
    res+="\t"+startStr;
    res+="\t"+endStr;
    res+="\t0";
    res+="\t1";
//cout<<"in get 3"<<endl;
    res+="\t"+to_string(static_cast<long long>(end-start));
    res+="\t0";
//cout<<"in get 4"<<endl;
    return res;
}

string Group::getAllMemBEDStr(Exons & es, ClusterGraph & cf, Clusters & ct, float p1, float p2, float p3)
{
//cout<<"in getAllMemBEDStr"<<endl;
    string res="";
    string mkey = members[0].getCkey();
    int numAlignM=members[0].getNumAlignments();
    vector<string> mAnnotStr=members[0].getAnnotStr(es);
    string mMultiCluster="";
    string mSplitCluster="";
    for(uint i=1;i<members.size();i++)
    {
//cout<<"i="<<i<<endl;
        int type = members[i].getMemberType(cf, es, mkey, numAlignM, mAnnotStr, p1, p2, p3);
//cout<<"type="<<type<<endl;
        if (type==1) //multi
        {   
            mMultiCluster=string_append(mMultiCluster,members[i].getCkey(),",");
        }
    }
    vector<string> multiStrVec;
    vector<string> splitStrVec;
    for(uint i=1;i<members.size();i++)
    {
//cout<<"i="<<i<<endl;
        int type = members[i].getMemberType(cf, es, mkey,numAlignM, mAnnotStr, p1, p2, p3);
//cout<<"type="<<type<<endl;
        if (type==1) //multi
        {
//cout<<"1111"<<endl;
            string multiStr = get_bed_12_str(members[i].getCkey(), ct);
            string relatedStr="M:"+mkey+","+mMultiCluster;
            multiStr+="\t"+to_string(static_cast<long long>(members[i].getNumAlignments()));
            multiStr+="\t"+relatedStr;
            multiStr+="\t"+members[i].getMemberStr2(es,cf,members[i].getCkey());
            multiStrVec.push_back(multiStr);
        }
        else if (type==2) //reads sharing split
        {
//cout<<"2222"<<endl;
            mSplitCluster=string_append(mSplitCluster,members[i].getCkey(),",");
            string splitStr = get_bed_12_str(members[i].getCkey(), ct);
            string relatedStr="F:"+mkey;
            splitStr+="\t"+to_string(static_cast<long long>(members[i].getNumAlignments()));
            splitStr+="\t"+relatedStr;
            splitStr+="\t"+members[i].getMemberStr2(es,cf,members[i].getCkey());
            splitStrVec.push_back(splitStr);
        }
    }
//cout<<"here 1"<<endl;
    string mAllRelatedStr="";
    if(mMultiCluster.compare("")!=0)
    {
        mAllRelatedStr="M:"+mkey+","+mMultiCluster;
    }
//cout<<"here 2"<<endl;
    if(mSplitCluster.compare("")!=0)
    {
        mSplitCluster="S:"+mSplitCluster;
        mAllRelatedStr=string_append(mAllRelatedStr,mSplitCluster,";");
    }
//cout<<"here 3"<<endl;
    if(mAllRelatedStr.compare("")==0)
        mAllRelatedStr="NA";
//cout<<"here 4"<<endl;
    string mainStr=get_bed_12_str(mkey,ct);
//cout<<"here 4.1"<<endl;
    mainStr+="\t"+to_string(static_cast<long long>(members[0].getNumAlignments()));
//cout<<"here 4.2"<<endl;
    mainStr+="\t"+mAllRelatedStr;
//cout<<"here 4.3"<<endl;
    mainStr+="\t"+getAllMemAbsorbStr(es,cf);
//cout<<"here 5"<<endl;
    res=res+mainStr;
    for(uint i=0;i<multiStrVec.size();i++)
    {
        res=res+"\n"+multiStrVec[i];
    }
//cout<<"here 6"<<endl;
    for(uint i=0;i<splitStrVec.size();i++)
    {
        res=res+"\n"+splitStrVec[i];
    }
//cout<<"here 7"<<endl;
    return res;    
}

void Groups::setGroups(vector<vector<string> > & strGroups, ClusterGraph & cg, Alignments & al, ReadAlignments & ra, Clusters & c, Exons & es)
{
    for(uint i=0;i<strGroups.size();i++)
    {
        Group g;
        g.setMembers(strGroups[i],cg,al,ra,c,es);
        g.sortMembers();//useful only for print
        groups.push_back(g);
    }

    c.moveToBegin();
    string ckey="";
    while(c.getNextKey(ckey))
    {
        if(cg.isClusterIn(ckey)==false)
        {
            vector<string> tmp;
            tmp.push_back(ckey);
            Group g;
            g.setMembers(tmp,cg,al, ra ,c,es);
            groups.push_back(g);
        }
    }

}

bool my_sort_group_func_1(Group i, Group j)
{
    int numi=i.getAllAlignmentSize();
    int numj=j.getAllAlignmentSize();
    if(i.getMembersSize()>j.getMembersSize())
        return true;
    else if(i.getMembersSize()==j.getMembersSize())
    {
        if(numi>numj)
            return true;
        else
            return false;
    }
    else
        return false;
}

void Groups::print(Exons & es, ClusterGraph & cf, string outputfile)
{
    ofstream of;
    of.open(outputfile);

    sort(groups.begin(),groups.end(),my_sort_group_func_1);
    for(uint i=0;i<groups.size();i++)
    {
        of<<endl;
        of<<"Group "<<i<<", member cluster number: ";
        of<<groups[i].getAllMemStr(es,cf);
    }
    of.close();
}

void Groups::print_bed(Exons & es, ClusterGraph & cf, Clusters & ct, string outputfile, string trackName, string description, float p1, float p2, float p3)
{
    ofstream of;
    of.open(outputfile);
    sort(groups.begin(),groups.end(),my_sort_group_func_1);
    of<<"track name=\""<<trackName<<"\" description=\""<<description<<"\" useScore=0";
    for(uint i=0;i<groups.size();i++)
    {
        of<<endl;
        of<<groups[i].getAllMemBEDStr(es,cf,ct,p1,p2,p3);
    }
    of.close();
}




