#include <string>


typedef struct
{
    string read_name;
    char strandType; //S for same, D for differnt, B for both
    string chr;
    int strand;
    uint32_t startPos;
    string sigar;
    int editDistance;
} alignment_t;


