#include "Utils.h"

std::vector<std::string> &my_split(const std::string &s, char delim, std::vector<std::string> &elems) {
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}

std::vector<std::string> my_split(const std::string &s, char delim) {
       std::vector<std::string> elems;
       my_split(s, delim, elems);
       return elems;
}

char targetChar;

bool is_same_char(char & currentChar)
{
    if (currentChar==targetChar)
        return true;
    else
        return false;
}

void removeChar(string & str, char a)
{
    targetChar=a;
    str.erase(remove_if(str.begin(), str.end(), is_same_char), str.end());
}



uint32_t getFilelength(char *file)
{
    struct stat filestatus;
    stat( file, &filestatus );
    uint32_t length = filestatus.st_size;
    return length;
}

int readBlock(char * block, int length, FILE *infile)
{
    bool a=fread(block,1,length,infile);
    if (a==false)
    {
        cerr<<"fread failed"<<endl;
        exit(1);
    }
    if(block[length-1]!='\n' && block[length-1]!=EOF)
    {
        bool b=fgets(block+length,1024,infile);
        if (b==false)
        {
            cerr<<"fread failed"<<endl;
            exit(1);
        }
        length+=strlen(block+length);
    }
    return length;
} 
