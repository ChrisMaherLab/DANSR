#include "Gene3.h"

void Exons::loadFromGTF(char * filename)
{
    uint32_t length=getFilelength(filename);
    FILE *infile;
    infile = fopen(filename, "r");
    if (!infile) {
        cerr<<"Couldn't open file for reading: "<<filename<<"."<<endl;
        exit(1);
    }

    char * fileContent;
    try
    {
        fileContent= new char [length];
    }
    catch(exception& e)
    {
        cerr << "Trying to allocate Memory to load Genes. exception caught: " << e.what() << endl;
        exit(1);
    }

    size_t result = fread (fileContent,1,length,infile);
    if (result != length)
    {
        cerr << "Fail to read the gene gtf file"<<endl;
        exit (1);
    }

    char* line_buffer=fileContent;
    char* nextN=NULL;
    unordered_set<string> tmpT;
    unordered_set<string> tmpG;

    while (1) {
        nextN=strchr(line_buffer,'\n');
        nextN[0]='\0';
        if(*line_buffer != '#')
        {
            vector<string> tmp_line = my_split(string(line_buffer), '\t');
            //int numnum=sscanf(line_buffer,"%s\t%s\t%s\t%u\t%u\t%s\t%s\t%s\t%s\n", chr, source, feature, &start, &end, tmp1, strand, tmp2, attribute);
            if(tmp_line.size()!=9)
            {
                cerr<<"Error loading genes at: "<<line_buffer<<endl;
                cerr<<"Not 9 columns."<<endl;
                exit(1);
            }
 
            if(tmp_line[2].compare("exon")==0)//feature
            {
                Exon ex;
                ex.setChr(tmp_line[0]);
                ex.setStart(atol(tmp_line[3].c_str()));
                ex.setEnd(atol(tmp_line[4].c_str()));
                if(tmp_line[6].compare("+")==0)
                    ex.setStrand(0);
                else if(tmp_line[6].compare("-")==0)
                    ex.setStrand(1);
                else
                    ex.setStrand(2);
                ex.setSource(tmp_line[1]);//source
                vector<string> tmp = my_split(tmp_line[8], ' ');//attribute
                for (uint i=0;i<tmp.size();i++)
                {
                    removeChar(tmp[i],'\"');
                    removeChar(tmp[i],';');
                }    
                int index_tId=-1;
                int index_gId=-1;
                int index_gName=-1;
                int index_gBioType=-1;
                for(uint i=0;i<tmp.size();i++)
                {
                    if(tmp[i].compare("transcript_id")==0)
                        index_tId=i+1;
                    if(tmp[i].compare("gene_id")==0)
                        index_gId=i+1;
                    if(tmp[i].compare("gene_name")==0)
                        index_gName=i+1;
                    if(tmp[i].compare("gene_biotype")==0)
                        index_gBioType=i+1;            
                }
                if (index_tId == -1 or index_gId==-1)
                {
                    cerr << "Attribute should contain transcript_id and gene_id."<<endl;
                    exit (1);
                }
                ex.setTId(tmp[index_tId]);
                ex.setGId(tmp[index_gId]);
                if(index_gName!=-1)
                    ex.setGName(tmp[index_gName]);
                else
                    ex.setGName(ex.getGId());
                if(index_gBioType!=-1)
                    ex.setGBioType(tmp[index_gBioType]);
                else
                    ex.setGBioType("NA");            

//cout<<ex.getTId();
//cout<<ex.getGId();
//cout<<ex.getGName();
//cout<<ex.getGBioType();
//cout<<endl;
		//ex.printExon(); //debugMJI
		if (1) //(ex.getGName().compare("Y_RNA"))
		{
		  exons.push_back(ex);
		  tmpT.insert(ex.getTId());
		  tmpG.insert(ex.getGId());
		}
            }
         }
         if(nextN-fileContent >= length-1)
         {
            cerr<<exons.size()<<" lines of exons loaded, "<<tmpT.size()<<" transcripts, "<<tmpG.size()<<" genes"<<endl;
            break;
         }
         else
         {
            line_buffer=nextN+1;
         }
    }
}

void Exons::loadGTFList(char * filenames)
{
    string tmpStr=string(filenames);
    removeChar(tmpStr, ' ');
    vector<string> tmp = my_split(tmpStr, ',');
    for(uint i=0;i<tmp.size();i++)
    {
      cout << "Processing GTF: " << tmp[i].c_str() << "\n"; //debugMJI
        loadFromGTF((char *)(tmp[i].c_str()));
    }
}

bool my_sort_exon_func (Exon i, Exon j)
{
    if(i.getChr().compare(j.getChr())<0)
    {
        return true;
    }
    else if(i.getChr().compare(j.getChr())==0)
    {
        if(i.getStart()<j.getStart())
            return true;
        else 
            return false;
    }
    else
    {
        return false;
    }
}

void Exons:: sortByCoodinate()
{
    sort(exons.begin(),exons.end(),my_sort_exon_func);
}

void Exons:: setMaps()
{
  cout << "exons.size() " << exons.size() << endl; //debugMJI
    for(uint i=0;i<exons.size();i++)
    {
      //exons[i].printExon(); //debugMJI
      /*
      vector<string> tKeys;
        geneToTranscripts(exons[i].getGId(),tKeys);
        bool exist=false;
	cout << "tKeys.size() " << tKeys.size() << endl; //debugMJI
        for(uint j=0;j<tKeys.size();j++)
        {
            if(tKeys[j].compare(exons[i].getTId())==0)
                exist=true;
        }
        if(exist==false)
            geneToTrans.insert(make_pair<string,string>(exons[i].getGId(),exons[i].getTId()));
      */
        transToExon.insert(make_pair(exons[i].getTId(),i));//keep all the exons
//cout<<exons[i].getTId()<<"\t"<<i<<endl;
    }
}


bool is_exon_overlapping(Exon & e1, Exon & e2)
{
    bool a = e1.getEnd() >= e2.getStart() && e2.getEnd() >= e1.getStart();
    bool b = e2.getEnd() >= e1.getStart() && e1.getEnd() >= e2.getStart();
    bool c = true;
    if (e1.getStrand()+e2.getStrand()==1)
        c=false;
    if ((a || b) && c)
        return true;
    return false;
}

void Exons:: coordinateToExons(string chr, uint32_t pos, uint32_t pos2, int strand, vector<uint> & ekeys)
{
    Exon dumbExon;
    dumbExon.setChr(chr);
    dumbExon.setStart(pos2);

    vector<Exon>::iterator up=upper_bound(exons.begin(),exons.end(),dumbExon,my_sort_exon_func);
    if (up-exons.begin()==0)
        return;
    up--;

    dumbExon.setStart(pos);
    dumbExon.setEnd(pos2);
    dumbExon.setStrand(strand);

    while(up-exons.begin()>=0 && (*up).getChr()==chr && pos2+10000 > (*up).getStart())
    {
        if(is_exon_overlapping(*up,dumbExon))
        {
            ekeys.push_back(up-exons.begin());
        }
        up--;
    }
}

/*
void Exons:: geneToTranscripts(string gkey, vector<string> & tkeys)
{
    auto range = geneToTrans.equal_range(gkey);
    unordered_multimap<string,string>::iterator it;
    for(it=range.first;it!=range.second;it++)
    {
        tkeys.push_back(it->second);
	//std::cout << "tkeys update; key: " << it->first <<" val: "<< it->second << "\n"; // debugMJI
    }
}
*/
void Exons:: transcriptToExons(string tkey, vector<uint> & eIds)
{
    auto range = transToExon.equal_range(tkey);
    unordered_multimap<string,uint>::iterator it;
    for(it=range.first;it!=range.second;it++)
    {
        eIds.push_back(it->second);
    }
}
