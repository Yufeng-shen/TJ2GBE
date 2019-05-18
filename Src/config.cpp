#include "config.hpp"

void Config::Tokenize(vector< vector<string> > & vsTokens, const string & sBuf, const string &sDelimiters){
    size_t iCurrentPos=0;
    size_t iDelimiterPos=0;
    size_t iEndofLine=0;
    while(iCurrentPos < sBuf.size()){
        vector<string> vsCurrentLine;
        iEndofLine = sBuf.find_first_of("\n",iEndofLine+1);
        if(iEndofLine == string::npos)
            return;
        while(iCurrentPos < iEndofLine){
            string currentToken;
            iDelimiterPos = sBuf.find_first_of(sDelimiters,iDelimiterPos);
            iCurrentPos= sBuf.find_first_not_of(sDelimiters,iCurrentPos);
            if((iCurrentPos == string::npos) || (iDelimiterPos == string::npos)) 
                break;
            for(size_t i= iCurrentPos;i<iDelimiterPos;i++)
                currentToken.append(sizeof(char),sBuf.at(i));
            if(currentToken.size()>0)
                vsCurrentLine.push_back(currentToken);
            iCurrentPos=iDelimiterPos+1;
            iDelimiterPos++;
        }
        vsTokens.push_back(vsCurrentLine);
    }
}

bool Config::Parse(const string & sBuf){
    const char* sKeywordStrings[]={
        "#", // comment
        "tripleJunctionFileName",
        "symmetryFileName",
        "numTJ",
        "maxNeighbor",
        "threshold",
        "fmax",
        "n",
        "ksym"
    };
    enum EKeyword{
        Ecomment,
        EtripleJunctionFileName,
        EsymmetryFileName,
        EnumTJ,
        EmaxNeighbor,
        Ethreshold,
        Efmax,
        En,
        Eksym,
        EnumKeyword // number of keywords, include comment
    };
    bool* bKeywordCheck;
    bKeywordCheck = new bool[EnumKeyword];
    for(int ii=0;ii<EnumKeyword;ii++)
        bKeywordCheck[ii]=false;
    bKeywordCheck[Ecomment]=true;
    vector< vector<string> > vsTokens;
    Tokenize(vsTokens, sBuf, " \t\n");
    for(size_t ii=0;ii<vsTokens.size();ii++){
        size_t iFirstToken;
        if(vsTokens[ii].size()==0)
            iFirstToken=Ecomment;
        else if (vsTokens[ii][0].find_first_of(sKeywordStrings[Ecomment])==0)
            iFirstToken=Ecomment;
        else{
            for(iFirstToken=0;iFirstToken<EnumKeyword;iFirstToken++){
                if(strcmp(sKeywordStrings[iFirstToken],vsTokens[ii][0].c_str())==0)
                    break;
            }
            switch(iFirstToken){
                case Ecomment:
                    break;
                case EtripleJunctionFileName:
                    tripleJunctionFileName = vsTokens[ii][1];
                    bKeywordCheck[EtripleJunctionFileName]=true;
                    break;
                case EsymmetryFileName:
                    symmetryFileName = vsTokens[ii][1];
                    bKeywordCheck[EsymmetryFileName]=true;
                    break;
                case EnumTJ:
                    numTJ = std::stoi(vsTokens[ii][1]);
                    bKeywordCheck[EnumTJ]=true;
                    break;
                case EmaxNeighbor:
                    maxNeighbor = std::stoi(vsTokens[ii][1]);
                    bKeywordCheck[EmaxNeighbor]=true;
                    break;
                case Ethreshold:
                    threshold = std::stod(vsTokens[ii][1]);
                    bKeywordCheck[Ethreshold]=true;
                    break;
                case Eksym:
                    ksym = std::stoi(vsTokens[ii][1]);
                    bKeywordCheck[Eksym]=true;
                    break;
                case Efmax:
                    for(int jj=0;jj<5;jj++)
                        fmax[jj]=std::stod(vsTokens[ii][jj+1]);
                    bKeywordCheck[Efmax]=true;
                    break;
                case En:
                    for(int jj=0;jj<5;jj++)
                        n[jj]=std::stoi(vsTokens[ii][jj+1]);
                    bKeywordCheck[En]=true;
                    break;
                default:{
                    std::cout<< "CfgFile Error: keyword not recognized: Line "<<ii<<" : "<<vsTokens[ii][0]<<std::endl;
                    return false;
                        }
            }
        }
    }
    bool notMissing=true;
    for(int ii=0;ii<EnumKeyword;ii++){
        if(!bKeywordCheck[ii]){
            std::cout<<"CfgFile Error: missing keyword "<<sKeywordStrings[ii]<<std::endl;
            notMissing=false;
        }
    }
    return notMissing;
}

bool Config::InputConfigParameters(const string & filename){
    std::ifstream inputFile(filename.c_str());
        if(inputFile.is_open()){
            std::stringstream buffer;
            buffer << inputFile.rdbuf();
            return Parse(buffer.str());
        }
        else{
            std::cout<<"Error opening: "<<filename<<std::endl;
            return false;
        }
}

void Config::PrintFile(){
    std::cout<<"triple junction filename: "<< tripleJunctionFileName<<std::endl
        <<"symmetry filename: "<<symmetryFileName<<std::endl
        <<"number of triple junctions: "<<numTJ<<std::endl
        <<"max number of neighbors: "<<maxNeighbor<<std::endl
        <<"threshold for neighbors: "<<threshold<<std::endl
        <<"ksym: "<<ksym<<std::endl
        <<"pi/(max Euler angle): "<<fmax[0]<<","<<fmax[1]<<","<<fmax[2]<<","<<fmax[3]<<","<<fmax[4]<<std::endl
        <<"number of divisions: "<<n[0]<<","<<n[1]<<","<<n[2]<<","<<n[3]<<","<<n[4]<<std::endl;
}
