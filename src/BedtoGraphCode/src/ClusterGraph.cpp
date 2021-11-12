/*
 * ClusterGraph.cpp
 *
 *  Created on: Mar 14, 2017
 *      Author: jinzhang
 */

#include "ClusterGraph.h"

using namespace std;

ClusterGraph::ClusterGraph() {
    // TODO Auto-generated constructor stub

}

ostream& operator<<( ostream &out,  ReadEdge &et )
{
    out<<et.getSize();
    return out;
}

ostream& operator<<( ostream &out, const ALGraph<int,ReadEdge> &g)
{

    int verNum = g.getVertexCount(),
        edgeNum = g.getEdgeCount();

    out << "This graph has " << verNum << " vertexes and " << edgeNum << " edges." << endl;
#if PRINTME
    for( int i=0; i<verNum; ++i )
    {
        int x1 = g.getData(i);
        out << x1 << " :    ";
        int j = g.getNextDst(x1);
        if( j != -1 )
        {
            int x2 = g.getData(j);
            out << "( " << x1 << ", " << x2 << ", " <<*(g.getWeight(x1,x2)) << " )" << "    ";
            do
            {
                j = g.getNextDst( x1, x2 );
                if( j != -1 )
                {
                    x2 = g.getData(j);
                    out << "( " << x1 << ", " << x2 << ", " << *(g.getWeight(x1,x2)) << " )" << "    ";
                }
                else
                    break;
            }
            while( j != -1 );
        }
        out << endl;
    }
#endif // PRINTME
    return out;
}



int ClusterGraph::printCg(Clusters & cs) {
//cs will be used to print detailed info for the nodes

    int verNum = fg.getVertexCount(),
        edgeNum = fg.getEdgeCount();

    cout << "This graph has " << verNum << " vertexes and " << edgeNum << " edges." << endl;

    return 0;
}

int ClusterGraph::printCg() {

    int verNum = fg.getVertexCount(),
        edgeNum = fg.getEdgeCount();

    cout << "This graph has " << verNum << " vertexes and " << edgeNum << " edges." << endl;

    return 0;
}


bool ClusterGraph::isClusterIn(string key) {
    return fg.vertexExists(key);
}

int ClusterGraph::getNeighbors(string key, vector<string> & ckeys){
    fg.neighboringVertexes(key, back_inserter(ckeys));
    return 0;
}

int ClusterGraph::addRead(string ckey1, string ckey2, string rKey) {

    ReadEdge *pt = fg.getWeight(ckey1,ckey2);
    ReadEdge *pt2 = fg.getWeight(ckey2,ckey1);

    if(pt==NULL && pt2==NULL)
    {
        ReadEdge re;
        re.addReadKey(rKey);
        fg.insertEdge(ckey1,ckey2,re);
        return 0;
    }
    else
    {
        pt->addReadKey(rKey);
        pt2->addReadKey(rKey);
        return 0;
    }

    return 0;
}


int ClusterGraph::updateReads(string ckey1, string ckey2, unordered_set<string> & rKeys) {

    ReadEdge * re = fg.getWeight(ckey1,ckey2);

    re->updateReadKeys(rKeys);

    return 0;
}

int ClusterGraph::removeRead(string ckey1, string ckey2, string rKey) {

    ReadEdge* re = fg.getWeight(ckey1,ckey2);
    if (re)
    {
        re->removeReadKey(rKey);    
        if (re->getSize()==0)
            fg.removeEdge(ckey1,ckey2);
    }
    return 0;
}

int ClusterGraph::removeEdge(string key1, string key2) {
    fg.removeEdge(key1,key2);
    return 0;
}

int ClusterGraph::getReadsNum(string key1, string key2) {
    int size=0;
    ReadEdge * re = fg.getWeight(key1,key2);
    if(re)
    {
        size=re->getSize();
    }
    return size;
}


int ReadEdge::addReadKey(string key)
{
    read_keys.insert(key);
    return 0;
}

int ReadEdge::removeReadKey(string key)
{
    unordered_set<string>::iterator it = read_keys.find(key);
    if (it!=read_keys.end())
    {
        read_keys.erase(it);
    } 
    return 0;
}

int ReadEdge::updateReadKeys(unordered_set<string> & rKeys)
{
    read_keys=rKeys;
    return 0;
}

int ReadEdge::getSize()
{
    return read_keys.size();
}

struct Connector {
    Connector(ClusterGraph & cg, vector<vector<string> > & res, unordered_set<string> & explored)
        : cg(cg)
        , explored(explored)
        , res(res)
    {
    }

    bool operator()(string ckey) {
//cout<<"ckey="<<ckey<<endl;
        if(explored.find(ckey)!=explored.end())
        {
//cout<<"in if"<<endl;
            return true;
        }
        else
        {
//cout<<"in else"<<endl;
            unordered_set<string> sub;
            vector<string> current;
            current.push_back(ckey);
            sub.insert(ckey);
            int size=1;           
            
            for(int ii=0;ii<size;ii++)
            {
//cout<<"loop"<<endl;
//cout<<"key="<<current[ii]<<endl;
                string key=current[ii];
                vector<string> tmp;
                cg.getNeighbors(key,tmp);
                for(uint x=0;x<tmp.size();x++)
                {
//cout<<"x="<<x<<endl;
                    string keyN=tmp[x];
                    if(explored.find(keyN)==explored.end())
                        explored.insert(keyN);
//cout<<"1111"<<endl;
                    if(sub.find(keyN)==sub.end())
                    {
//cout<<"in this if"<<endl;
                        current.push_back(keyN);
                        size++;
                        sub.insert(keyN);
//cout<<"2222"<<endl;
                    }
                }
            }
//cout<<"before push"<<endl;
            res.push_back(current);
            return true;
        }
    }   
    
    ClusterGraph& cg;
    unordered_set<string> & explored;
    vector<vector<string> > & res;

};


struct Connector2 {
    Connector2(ClusterGraph & cg, Clusters & cs, vector<vector<string> > & res, unordered_set<string> & in_sub_or_self, float percent_cutoff)
        : cg(cg)
        , cs(cs)
        , in_sub_or_self(in_sub_or_self)
        , res(res)
        , percent_cutoff(percent_cutoff)
    {
    }

    bool operator()(string ckey) {
        if(in_sub_or_self.find(ckey)!=in_sub_or_self.end())
        {
            return true;
        }
        else
        {
            unordered_set<string> sub;
            vector<string> current;
            current.push_back(ckey);
            sub.insert(ckey);
            int size=1;

            for(int ii=0;ii<size;ii++) // while might be clearer
            {
                string key=current[ii];
                vector<string> tmp;
                cg.getNeighbors(key,tmp);
                for(uint x=0;x<tmp.size();x++)
                {
                    string keyN=tmp[x];
                    //code for checking whether pass cut off
                    int numReads=cg.getReadsNum(key,keyN);
                    bool is_pass = false;
                    int numA = cs.numAlignments(key);
                    int numB = cs.numAlignments(keyN);
                    if (numReads > numA*percent_cutoff || numReads > numB*percent_cutoff)
                        is_pass=true; 
                    if(is_pass==true)
                    {
                        if(in_sub_or_self.find(key)==in_sub_or_self.end())
                            in_sub_or_self.insert(key);
                        if(in_sub_or_self.find(keyN)==in_sub_or_self.end())
                            in_sub_or_self.insert(keyN);
                        if(sub.find(keyN)==sub.end())
                        {
                            current.push_back(keyN);
                            size++;
                            sub.insert(keyN);
                        }
                    }
                }
            }
            if(size==1)
                in_sub_or_self.insert(ckey);
            res.push_back(current);
            return true;
        }
    }

    ClusterGraph& cg;
    Clusters& cs;
    unordered_set<string> & in_sub_or_self;
    vector<vector<string> > & res;
    float percent_cutoff;
};


int ClusterGraph::getDegree(string ckey)
{
    vector<string> nei;
    getNeighbors(ckey, nei);
    return nei.size();    
}


void ClusterGraph::getConnectedClusters(vector<vector<string> > & res)
{
    unordered_set<string> explored;
    Connector connector(*this, res, explored);
    fg.foreachVertex(connector);
}

void ClusterGraph::getConnectedClusters(Clusters & cs, vector<vector<string> > & res, float percent_cutoff)
{
    unordered_set<string> in_sub_or_self;
    Connector2 connector2(*this, cs, res, in_sub_or_self, percent_cutoff);
    fg.foreachVertex(connector2);
}



