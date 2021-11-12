#include "BedToGraph.h"

void BedToGraph::connectBedFile(char * filename)
{
    bedfile.open (filename, std::ifstream::in);
}

bool BedToGraph::getNextBedLine(string & line)
{
    getline(bedfile,line);
    return (bool)bedfile;
}

void BedToGraph::closeConnection()
{
    bedfile.close();
}

void BedToGraph::loadAlignmentsAndClusters(Alignments & a, Clusters & c, AlignmentClusters & ac, ReadAlignments & ra)
{
//cout<<"in loadAlignmentsAndClusters"<<endl;
    string line;
    while(getNextBedLine(line))
    {
//cout<<"line:"<<line<<endl;
        Cluster cc;
        if(cc.changeFromLine(line,bedFlag)) //write this
        {
//cout<<"if 1"<<endl;
            string cKey=cc.getKey();
            vector<string> als;
            if(a.getAlignmentStrings(line,als))  //write this
            {
//cout<<"if 2"<<endl;
                for(uint x=0;x<als.size();x++)
                {
//cout<<"for 3"<<endl;
                    Alignment aa;
                    if(aa.changeFromString(als[x])) //write this
                    {
//cout<<"if 4"<<endl;
                        string aKey=aa.getKey();
                        a.insertAlignment(aKey,aa);
                        ac.insert(aKey,cKey);
//cout<<"insert "<<aKey<<" "<<cKey<<endl;
                        ra.insert(aa.getReadName(),aKey);
                        cc.addAlignmentKey(aKey);
//cout<<"insert "<<aa.getReadName()<<" "<<aKey<<endl;
                    }
                }
            }
            c.insertCluster(cKey,cc);
        }
    }
    closeConnection();
}

void BedToGraph::createGraph(ClusterGraph & cg, ReadAlignments & ra, Alignments & a, Clusters & c, AlignmentClusters & ac)
{
    ra.moveToBegin();
    string readName="";
    while(ra.getNextName(readName))
    {
//cout<<"Name is: "<<readName<<endl;
        vector<string> aKeys;
        ra.getAlignmentKeys(readName,aKeys);
        vector<string>::iterator it;
        vector<string> cKeys;
        for(it=aKeys.begin();it!=aKeys.end();it++)
        {
            vector<string> tmp;
            string akey=(*it);
            ac.getClusterKeys(akey,tmp);
            cKeys.insert(cKeys.end(),tmp.begin(),tmp.end());
        }
        for(uint x=0;x<cKeys.size();x++)
        {
            for(uint y=x+1;y<cKeys.size();y++)
            {
//cout<<"add "<<cKeys[x]<<" "<<cKeys[y]<<" "<<readName<<endl;
                 if (cKeys[x].compare(cKeys[y])!=0)    // diffrent alignments of the same read can be in the same cluster
                     cg.addRead(cKeys[x], cKeys[y],readName);
            }
        }
    }
}



