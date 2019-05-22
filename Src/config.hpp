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
        string tripleJunctionFileName; // filename of the triple junctions
        string symmetryFileName; // filename of the symmetry file
        string outputDir; // output directory for saving results
        unsigned int numTJ; // number of triple junctions
        int maxNeighbor; // number of maximum similar grain boundaries
        double threshold; // threshold of being a similar grain boundary
        double fmax[5]; // subdomain size
        int n[5]; // number of divisions of each dimension of subdomain
        int ksym; 

        bool InputConfigParameters(const string & filename);
        void PrintFile();
};


#endif
