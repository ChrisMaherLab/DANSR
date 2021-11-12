/*
 * Clusters.h
 *
 *  Created on: March 15, 2017
 *      Author: jinzhang
 *
 */

#ifndef CLUSTERS_H_
#define CLUSTERS_H_

#include <iostream>
#include <unordered_set>
#include <unordered_map>
#include <vector>
#include <iterator>
#include <cstdlib>
#include <cstdint>

using namespace std;

#include "Utils.h"
#include "Alignments.h"
#include "ReadAlignments.h"

class Cluster
{
private:
    unordered_set<string> alignment_keys;
    string chr;
    uint32_t start;
    int length;
    int strand;
    int uniq;   
 
public:
    Cluster(){};
    void setChr(string chr){this->chr=chr;}
    string getChr(){return this->chr;}
    void setStart(uint32_t start){this->start=start;}
    uint32_t getStart(){return this->start;}
    void setLength(uint32_t length){this->length=length;}
    int getLength(){return this->length;}
    void setStrand(int strand){this->strand=strand;}
    int getStrand(){return this->strand;}
    
    int numUniq(Alignments & as, ReadAlignments & ra);

    void addAlignmentKey(string key);
    bool isAlignmentIn(string key);
    void removeAlignmentKey(string key);
    string getKey();
    bool changeFromLine(string line, int bedFlag);//bedFlag is strand
    int numAlign(){return alignment_keys.size();}
};

class Clusters
{
private:
    unordered_map<string,Cluster> clusters;
    unordered_map<string,Cluster>::iterator cit;
public:
    Clusters(){};
    Cluster getCluster(string key);
    void insertCluster(string key, Cluster c);
    bool isClusterIn(string key);
    void moveToBegin(){this->cit=clusters.begin();}
    bool getNextKey(string & ckey);
    int size(){return clusters.size();}
    int numAlignments(string ckey);
    int numUniqAlignments(string ckey, Alignments & as, ReadAlignments & ra);
};

#endif
