#include <iostream>

using namespace std;

#include "Clusters.h"
#include "Alignments.h"
#include "ClusterGraph.h"
#include "Clusters.h"
#include "ReadAlignments.h"
#include "AlignmentClusters.h"
#include "BedToGraph.h"
#include "Gene3.h"
#include "Groups.h"

void getGeneIds(string ckey, Clusters & c, Exons & es, vector<string> & gIds)
{
  //cout<<"in here"<<endl;
  //cout<<ckey<<"######"<<endl;
    Cluster a=c.getCluster(ckey);
    vector<uint> ekeys;
    //cout<<"in here 1111"<<endl;
    es.coordinateToExons(a.getChr(),a.getStart(),a.getStart()+a.getLength()-1,a.getStrand(),ekeys);
    unordered_set<string> ids;
    //cout<<"in here 2222"<<endl;
    for(uint i=0;i<ekeys.size();i++)
    {
      //cout<<"insert "<<es.getExon(ekeys[i]).getGId()<<endl;
        ids.insert(es.getExon(ekeys[i]).getGId());
    }
    //cout<<"set size="<<ids.size()<<endl;
    unordered_set<string>::iterator it;
    for(it=ids.begin();it!=ids.end();it++)
        gIds.push_back(*it);
    //cout<<"in here 4444"<<endl;
}

/*
string getGstr(vector<string> & gIds, Exons & es)
{
  //cout<<"in here 2"<<endl;
    string res;
    for(uint i=0;i<gIds.size();i++)
    {
      //cout<<"gid="<<gIds[i]<<endl;
        vector<string> tkeys;
        es.geneToTranscripts(gIds[i], tkeys);
      //cout<<"tkeys.size="<<tkeys.size()<<endl;
      //for(uint x=0;x<tkeys.size();x++)
      //cout<<tkeys[x]<<endl;
      //cout<<"here XXXX"<<endl;
      //cout<<"eIds.size="<<eIds.size()<<endl;
        Exon e0;
        vector<uint32_t> left;
        vector<uint32_t> right;
        for(uint j=0;j<tkeys.size();j++)
        {
	  //cout<<endl;
	  //cout<<"tid="<<tkeys[j]<<endl;
            vector<uint> eIds;
            es.transcriptToExons(tkeys[j],eIds);
            for(uint k=0;k<eIds.size();k++)
            {
                e0=es.getExon(eIds[k]);
                left.push_back(e0.getStart());
                right.push_back(e0.getEnd());
		//cout<<endl;
		//cout<<e0.getGName()<<endl;
		//cout<<e0.getTId()<<endl;
		//cout<<e0.getGId()<<endl;
            }
        }
      //cout<<"here YYY"<<endl;
      //        Exon e0=es.getExon(eIds[0]);
      //cout<<"here ZZZ"<<endl;
      //cout<<e0.getGName()<<endl;
      //cout<<e0.getGBioType()<<endl;
      //cout<<e0.getChr()<<endl;
      //cout<<to_string(static_cast<long long>(*min_element(left.begin(),left.end())))<<endl;
      //cout<<to_string(static_cast<long long>(*max_element(right.begin(),right.end())))<<endl;
    
        res=res+e0.getGName()+"_"+e0.getGBioType()+"_"+e0.getChr()+"_"+to_string(static_cast<long long>(*min_element(left.begin(),left.end())))+"_"+to_string(static_cast<long long>(*max_element(right.begin(),right.end())))+"|";
    }
    //cout<<"here 000"<<endl;
    //cout<<"left here 2"<<endl;
    //cout<<"res="<<res<<endl;
    return res;
}
*/

int main(int arg, char * argv[])
{
  
cout<<"in main"<<endl;
 
    Alignments a;
    Clusters c;
    ClusterGraph cg;
    ReadAlignments ra;
    AlignmentClusters ac;    
    BedToGraph btg;

    cout<<"load bed"<<endl;
     
    btg.connectBedFile(argv[1]);
    btg.setBedFlag(0);
    btg.loadAlignmentsAndClusters(a,c,ac,ra);
    btg.connectBedFile(argv[2]);
    btg.setBedFlag(1);
    btg.loadAlignmentsAndClusters(a,c,ac,ra);
    btg.createGraph(cg,ra,a,c,ac);


    //Test Gene3
    cout<<"load gtf"<<endl;
    Exons es;
    es.loadGTFList(argv[3]);
    es.sortByCoodinate();
    es.setMaps();

/*
    vector<uint> ekeys;
    //21	havana	exon	38474027	38474121
    //ENSG00000157554"; gene_version "18"; transcript_id "ENST00000481609"
    cout<<"coordinateToExons"<<endl;
    es.coordinateToExons("11",7110166,7110180,0,ekeys);
    for(uint i=0;i<ekeys.size();i++)
        cout<<ekeys[i]<<endl;
    cout<<endl;
    vector<string> tkeys;
    cout<<"geneToTranscripts"<<endl;
    es.geneToTranscripts("ENSG00000170748", tkeys);
    for(uint i=0;i<tkeys.size();i++)
        cout<<tkeys[i]<<endl;
    cout<<endl;
    ekeys.clear();
    cout<<"transcriptToExons"<<endl;
    es.transcriptToExons("ENST00000306904", ekeys);
    for(uint i=0;i<ekeys.size();i++)
        cout<<ekeys[i]<<endl;
    cout<<endl;
*/
    cg.printCg();

    cout<<"cluster size is: "<<c.size()<<endl;

    vector<vector<string> > res;
    //cg.getConnectedClusters(res);
    cg.getConnectedClusters(c, res, 0.5);
/*
    for(uint i=0;i<res.size();i++)
    {
        cout<<res[i].size()<<"\t";
        for(uint j=0;j<res[i].size();j++)
        {
//cout<<"res[i][j]:"<<res[i][j]<<endl;
            vector<string> gIds;
            getGeneIds(res[i][j],c,es,gIds);
//cout<<"gIds.size()="<<gIds.size()<<endl;
            cout<<res[i][j]<<"("<<cg.getDegree(res[i][j])<<")"<<"<"<<c.numAlignments(res[i][j])<<">"<<"["<<getGstr(gIds,es)<<"]; ";
        }
        cout<<endl;
    }

    cout<<"for clusters not connected"<<endl;
    c.moveToBegin();
    string ckey="";
    while(c.getNextKey(ckey))
    {
        if(cg.isClusterIn(ckey)==false)
        {
            vector<string> gIds;
            getGeneIds(ckey,c,es,gIds);
            cout<<"1\t"<<ckey<<"("<<0<<")"<<"<"<<c.numAlignments(ckey)<<">"<<"["<<getGstr(gIds,es)<<"]"<<endl;
        }
    }
*/    
    Groups gps;
    gps.setGroups(res, cg, a, ra, c, es);

    gps.print(es,cg,string(argv[4]));
    gps.print_bed(es,cg,c,string(argv[5]),"trackName", "description", 0.5, 0.1, 0.05);


/*
    //test container uniqueness:
    unordered_set<int> useta;
    useta.insert(1);useta.insert(1);useta.insert(1);
    cout<<useta.size()<<endl;
    unordered_map<int,int> usetb;
    usetb.insert(make_pair<int,int>(1,1));usetb.insert(make_pair<int,int>(1,1));usetb.insert(make_pair<int,int>(1,1));
    cout<<usetb.size()<<endl;
    unordered_multimap<int,int> usetc;
    usetc.insert(make_pair<int,int>(1,1));usetc.insert(make_pair<int,int>(1,1));usetc.insert(make_pair<int,int>(1,1));
    cout<<usetc.size()<<endl;

*/
    
    return 0;
}
