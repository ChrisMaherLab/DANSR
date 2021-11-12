/*
 * AlignmentClusters.h
 *
 *  Created on: March 16, 2017
 *      Author: jinzhang
 *
 */

#ifndef ALIGNMENTCLUSTERS_H_
#define ALIGNMENTCLUSTERS_H_

#include <iostream>
#include <unordered_set>
#include <unordered_map>
#include <vector>
#include <iterator>
#include <cstdint>

using namespace std;

class AlignmentClusters
{
private:
    unordered_multimap<string,string> alignmentClusters;
public:
    AlignmentClusters(){};
    void getClusterKeys(string akey, vector<string> & keys);
    void insert(string akey, string ckey);
};

#endif
