#ifndef TJ2GBE_CONFIG_HPP
#define TJ2GBE_CONFIG_HPP

#include <string>
#include <cstring>
#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>

using std::string;
using std::vector;

class Config{
    private:
        bool Parse(const string & sBuf);
        void Tokenize(vector< vector<string> > & vsTokens, const string & sBuf, const string & sDelimiters); // read the whole file sBuf to lines of tokens in vsTokens
    public:
        string tripleJunctionFileName;
        string symmetryFileName;
        unsigned int numTJ;
        int maxNeighbor;
        double threshold;
        double fmax[5];
        int n[5];
        int ksym;

        bool InputConfigParameters(const string & filename);
        void PrintFile();
};


#endif
