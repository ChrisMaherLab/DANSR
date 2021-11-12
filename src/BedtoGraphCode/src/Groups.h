/*
 * Groups.h
 *
 *  Created on: Apr 5 2017
 *      Author: jinzhang
 *
 */

#ifndef GROUPS_H_
#define GROUPS_H_

#include <iostream>
#include <list>
#include <unordered_set>
#include <vector>
#include <string>
#include <iterator>
#include <algorithm>
#include <fstream>

using namespace std;

#include "Clusters.h"
#include "Alignments.h"
#include "ClusterGraph.h"
#include "Clusters.h"
#include "ReadAlignments.h"
#include "AlignmentClusters.h"
#include "BedToGraph.h"
#include "Gene3.h"
#include "Utils.h"

class Member
{
private:
    string ckey;
    int degree;
    int numAlignments;
    int numUniqAlignments;
    vector<uint> eIds;
    int status;// 0 initial 1 absorbed
public:
    Member(){}
    string getCkey(){return this->ckey;}
    void setCkey(string ckey){this->ckey=ckey;}
    int getDegree(){return this->degree;}
    void setDegree(int degree){this->degree=degree;}
    void setDegree(ClusterGraph & cg){this->degree=cg.getDegree(ckey);}
    int getNumAlignments(){return this->numAlignments;}
    void setNumAlignments(int numAlignments){this->numAlignments=numAlignments;}
    void setNumAlignments(Clusters & c){this->numAlignments=c.numAlignments(ckey);}
    int getNumUniqAlignments(){return this->numUniqAlignments;}
    void setNumUniqAlignments(int numUniqAlignments){this->numUniqAlignments=numUniqAlignments;}
    void setNumUniqAlignments(Clusters & c, Alignments & al, ReadAlignments & ra){this->numUniqAlignments=c.numUniqAlignments(ckey,al, ra);}
    int getStatus(){return this->status;}
    void setStatus(int status){this->status=status;}
    void setEIds(vector<uint> eIds){this->eIds=eIds;}
    vector<uint> getEIds(){return this->eIds;}
    string getMemberStr(Exons & es, ClusterGraph & cf, string mkey);
    string getMemberStr2(Exons & es, ClusterGraph & cf, string mkey);
    vector<string> getAnnotStr(Exons & es);
    int getMemberType(ClusterGraph & cf, Exons & es, string mkey, int numAlignM, vector<string> mAnnotStr, float p1, float p2, float p3);

};

class Group
{
private:
    vector<Member> members;
public:
    Group(){};
    void setMembers(vector<string> & ckeys, ClusterGraph & cg, Alignments & al, ReadAlignments & ra ,Clusters & c, Exons & es);
    void sortMembers();
    int getMembersSize(){return members.size();}
    int getAllAlignmentSize();
    string getAllMemStr(Exons & es, ClusterGraph & cf);
    string getAllMemAbsorbStr(Exons & es, ClusterGraph & cf);
    string getAllMemBEDStr(Exons & es, ClusterGraph & cf, Clusters & ct, float p1, float p2, float p3);
};

class Groups
{
private:
    vector<Group> groups;
public:
    Groups(){};
    void setGroups(vector<vector<string> > & strGroups, ClusterGraph & cg, Alignments & al, ReadAlignments & ra, Clusters & c, Exons & es);
    void print(Exons & es, ClusterGraph & cf, string outputfile);
    void print_bed(Exons & es, ClusterGraph & cf, Clusters & ct, string outputfile,string trackName, string description,float p1,float p2, float p3);
};


#endif
