/*
 * Utils.h
 *
 *  Copied from old code:
 *      on: Mar 17, 2017
 *      Author: jinzhang
 */


#ifndef UTILS_H_
#define UTILS_H_

#include <iostream>
#include <sstream>
#include <string>
#include <cstring>
#include <vector>
#include <algorithm>
#include <sys/types.h>
#include <sys/stat.h>
#include <stdint.h>

using namespace std;

std::vector<std::string> &my_split(const std::string &s, char delim, std::vector<std::string> &elems);

std::vector<std::string> my_split(const std::string &s, char delim);

void removeChar(string & str, char a);

uint32_t getFilelength(char *file);
int readBlock(char * block, int length, FILE *infile);

#endif
