/*
 *  BedToGraph.h
 *
 *  Created on: March 16, 2017
 *      Author: jinzhang
 *
 */

#ifndef BEDTOGRAPH_H_
#define BEDTOGRAPH_H_

#include <iostream>
#include <iterator>
#include <fstream>
#include <string>

using namespace std;

#include "Alignments.h"
#include "ClusterGraph.h"
#include "Clusters.h"
#include "ReadAlignments.h"
#include "AlignmentClusters.h"

class BedToGraph
{
private:
    string currentBedLine;
    ifstream bedfile;
    int bedFlag;
public:
    void setBedFlag(int bedFlag){this->bedFlag=bedFlag;}
    int getBedFlag(){return this->bedFlag;}
    void connectBedFile(char * filename);
    bool getNextBedLine(string & line);
    void closeConnection();
    // add alignment,clusters,AlignmentClusters, ReadAlignment
    void loadAlignmentsAndClusters(Alignments & a, Clusters & c, AlignmentClusters & ac, ReadAlignments & ra);
    void createGraph(ClusterGraph & cg, ReadAlignments & ra, Alignments & a, Clusters & c, AlignmentClusters & ac); //Graph
};

#endif
