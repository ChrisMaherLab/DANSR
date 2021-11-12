/*
 * ReadAlignments.h
 *
 *  Created on: March 16, 2017
 *      Author: jinzhang
 *
 */

#ifndef READALIGNMENTS_H_
#define READALIGNMENTS_H_

#include <iostream>
#include <unordered_set>
#include <unordered_map>
#include <vector>
#include <iterator>
#include <algorithm>
#include <cstdint>

using namespace std;

class ReadAlignments
{
private:
    unordered_multimap<string,string> readAlignments;
    unordered_set<string> readNames;
    unordered_set<string>::iterator rit;
public:
    ReadAlignments(){};
    void moveToBegin(){rit=readNames.begin();}
    void getAlignmentKeys(string readName, vector<string> & keys);
    void insert(string readName, string alignmentKey);
    bool getNextName(string & readName);
};

#endif
