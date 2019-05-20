#include "subdomain.hpp"

void Subdomain::initialize(const double* fmax_, const int* n_,
        unsigned int numGB,  int numSym){
    for(int ii=0;ii<5;ii++){
        fmax[ii]=fmax_[ii];
        ndiv[ii]=n_[ii];
        afrange[ii]=fmax[ii];
    }
    afrange[1]=1-cos(afrange[1]);
    afrange[3]=1-cos(afrange[3]);
    ksym=numSym;
    this -> numGB = numGB;
    totaln=ndiv[0]*ndiv[1]*ndiv[2]*ndiv[3]*ndiv[4];
}

unsigned int Subdomain::neuler_to_cell(const double* af){
    int m[5];
    for(int ii=0;ii<5;ii++){
        m[ii]=std::min(int(af[ii]*ndiv[ii]/afrange[ii]),ndiv[ii]-1);
    }
    unsigned int res=m[0];
    for(int ii=1;ii<5;ii++){
        unsigned int tmp=m[ii];
        for(int jj=0;jj<ii;jj++){
            tmp*=ndiv[jj];
        }
        res+=tmp;
    }
    return res;
}

bool Subdomain::in_subdomain(const double* af){
    bool testor;
    testor=(af[0]<=fmax[0])&&(af[1]<=fmax[1])&&(af[2]<=fmax[2])&&(af[3]<=fmax[3])&&(af[4]<=fmax[4]);
    return testor;
}


void Subdomain::count_each_cell(const unsigned int* idxcell,
        const int* Qs){
    num = new unsigned int [totaln]();
    cumnum = new unsigned int [totaln];
    unsigned int s=ksym*numGB;
    for(unsigned int ii=0;ii<s;ii++){
        num[idxcell[ii]]+=1;
    }
    cumnum[0]=0;
    for(unsigned int ii=1;ii<totaln;ii++){
        cumnum[ii]=cumnum[ii-1]+num[ii-1];
    }
    int* offset = new int [totaln]();
    ids = new unsigned int [numGB*ksym];
    qs = new int [numGB*ksym];
    for(unsigned int ii=0;ii<numGB;ii++){
        for(int jj=0;jj<ksym;jj++){
            unsigned int idx= idxcell[ii*ksym+jj];
            ids[cumnum[idx]+offset[idx]]=ii;
            qs[cumnum[idx]+offset[idx]]=Qs[ii*ksym+jj];
            offset[idx]+=1;
        }
    }
    delete offset;
}

Subdomain::~Subdomain(){
    delete num;
    delete cumnum;
    delete ids;
    delete qs;
}

