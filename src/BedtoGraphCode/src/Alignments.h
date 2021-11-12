/*
 * Alignments.h
 *
 *  Created on: March 16, 2017
 *      Author: jinzhang
 *
 */

#ifndef ALIGNMENTS_H_
#define ALIGNMENTS_H_

#include <iostream>
#include <unordered_set>
#include <unordered_map>
#include <vector>
#include <iterator>
#include <cstdlib>
#include <cstdint>

using namespace std;

#include "Utils.h"

class Alignment
{
private:
    string readName;
    char strandType; //S for same, D for differnt, B for both
    string chr;
    int strand;
    uint32_t start;
    string cigar;
    int editDistance;

public:
    Alignment(){};
    void setReadName(string readName){this->readName=readName;}
    string getReadName(){return this->readName;}
    void setStrandType(char strandType){this->strandType=strandType;}
    char getStrandType(){return this->strandType;}
    void setChr(string chr){this->chr=chr;}
    string getChr(){return this->chr;}
    void setStrand(int strand){this->strand=strand;}
    int getstrand(){return this->strand;}
    void setStart(uint32_t start){this->start=start;}
    uint32_t getStart(){return this->start;}
    void setCigar(string cigar){this->cigar=cigar;}
    string getCigar(){return this->cigar;}
    void setEditDistance(int strand){this->strand=editDistance;}
    int getEditDistance(){return this->editDistance;}

    string getKey();

    bool changeFromString(string str);

};

class Alignments
{
private:
    unordered_map<string,Alignment> alignments;

public:
    Alignments(){};
    Alignment getAlignment(string key);
    void insertAlignment(string key, Alignment a);
    bool isAlignmentIn(string key);
    bool getAlignmentStrings(string line, vector<string> & strs);
};

#endif



