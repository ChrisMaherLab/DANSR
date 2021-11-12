/*
 * ClusterGraph.h
 *
 *  Created on: March 14, 2017
 *      Author: jinzhang
 *  
 *  Using the template of ALGraph, 
 */

#ifndef CLUSTERGRAPH_H_
#define CLUSTERGRAPH_H_

#include <iostream>
#include <list>
#include <unordered_set>
#include <vector>
#include <string>
#include <iterator>
#include <algorithm>

using namespace std;

#include "ALGraph.h"
#include "MyTypes.h"
#include "Clusters.h"

class ReadEdge{
    private:
        unordered_set<string> read_keys;
    public:
        int addReadKey(string key);
        int removeReadKey(string key);
        int updateReadKeys(unordered_set<string> & rKeys);
        int getSize();
};



class ClusterGraph {
private:

    typedef ALGraph<string, ReadEdge> GraphType;
    typedef GraphType::EdgeList EdgeList;
    typedef GraphType::iterator vertex_iterator;
    typedef EdgeList::iterator edge_iterator;
    typedef GraphType::const_iterator const_vertex_iterator;
    typedef EdgeList::const_iterator const_edge_iterator;
    
    GraphType fg;

public:
    
    ClusterGraph();
    int printCg(Clusters & cs);
    int printCg();
    bool isClusterIn(string key);

    int addRead(string ckey1, string ckey2, string rKey); 
    int removeRead(string ckey1, string ckey2, string rKey);
    
    int getNeighbors(string key, vector<string> & ckeys);

    int updateReads(string ckey1, string ckey2, unordered_set<string> & rkeys);
    int removeEdge(string key1, string key2);

    int getReadsNum(string key1, string key2);

    int getDegree(string ckey);

    void getConnectedClusters(vector<vector<string> > & res); 
    void getConnectedClusters(Clusters & cs, vector<vector<string> > & res, float percent_cutoff);

};

#endif /* CLUSTERGRAPH_H_ */
