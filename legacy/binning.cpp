#include "binning.hpp"
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>


inline void MATMAT_d4(const double* A, const double* B,
        double* res){
    for(int ii=0;ii<4;ii++){
        for(int jj=0;jj<4;jj++){
            res[4*ii+jj]=0;
            for(int kk=0;kk<4;kk++){
                res[4*ii+jj]+=A[4*ii+kk]*B[4*kk+jj];
            }
        }
    }
}
inline void MATVEC_d3(const double* A, const double* V,
        double* res){
    for(int ii=0;ii<3;ii++){
        res[ii]=0;
        for(int jj=0;jj<3;jj++){
            res[ii]+=A[3*ii+jj]*V[jj];
        }
    }
}
inline void Euler2Mat(const double* f, double* g){
    double sf=sin(f[1]);
    double cf=cos(f[1]);
    double sf1=sin(f[0]);
    double cf1=cos(f[0]);
    double sf2=sin(f[2]);
    double cf2=cos(f[2]);
    g[0*3+0]=cf1*cf2-sf1*sf2*cf;
  
    g[0*3+1]=sf1*cf2+cf1*sf2*cf;
  
    g[0*3+2]=sf2*sf;
   
    g[1*3+0]=-cf1*sf2-sf1*cf2*cf;
  
    g[1*3+1]=-sf1*sf2+cf1*cf2*cf;
  
    g[1*3+2]=cf2*sf;
   
    g[2*3+0]=sf1*sf;
  
    g[2*3+1]=-cf1*sf;
  
    g[2*3+2]=cf;
}
inline void TRANSPOSE_d3(const double* A, double* res){
    for(int ii=0;ii<3;ii++){
        for(int jj=0;jj<3;jj++){
            res[3*ii+jj]=A[3*jj+ii];
        }
    }
}
inline void TRANSPOSE_d4(const double* A, double* res){
    for(int ii=0;ii<4;ii++){
        for(int jj=0;jj<4;jj++){
            res[4*ii+jj]=A[4*jj+ii];
        }
    }
}

inline double acoss(double ang){
    if(ang>1) return 0.0;
    if(ang<-1) return PI;
    return acos(ang);
}

inline double distance(const double* A, const double* B){
    double res=0;
    for(int ii=0;ii<16;ii++){
        res += (A[ii]-B[ii])*(A[ii]-B[ii]);
    }
    return res;
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

Subdomain::~Subdomain(){
//    if(num!=NULL){
//    delete num;
//    delete cumnum;
//    delete ids;
//    delete qs;}
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

void TJ::normal_transformer(int q, const double* b, double* M){
    q=q-1;
    bool isignum=false;
    bool istranspose=false;
    if(q>=2*nselau*nselau){
        q=q-2*nselau*nselau;
        istranspose=true;
    }
    if(q>=nselau*nselau){
        q=q-nselau*nselau;
        isignum=true;
    }
    int i,j;
    i=int(q/nselau);
    j=q % nselau;

    if((!isignum)&&(!istranspose)){
        for(int ii=0;ii<3;ii++){
            for(int jj=0;jj<3;jj++){
                M[ii*3+jj]=gcsym[i*9+ii*3+jj];
            }
        }
    }
    else if((isignum)&&(!istranspose)){
        for(int ii=0;ii<3;ii++){
            for(int jj=0;jj<3;jj++){
                M[ii*3+jj]=-gcsym[i*9+ii*3+jj];
            }
        }
    }
    else if((!isignum)&&(istranspose)){
        double tmp[9];
        for(int iii=0;iii<3;iii++){
            for(int jjj=0;jjj<3;jjj++){
                tmp[iii*3+jjj]=-b[jjj*4+iii];
            }
        }
        MATMAT_d3(&gcsym[i*9],tmp,M);
    }
    else{
        double tmp[9];
        for(int iii=0;iii<3;iii++){
            for(int jjj=0;jjj<3;jjj++){
                tmp[iii*3+jjj]=b[jjj*4+iii];
            }
        }
        MATMAT_d3(&gcsym[i*9],tmp,M);
    }
}

void TJ::q2bx(int q, const double* b, double* bx){
    if(q<1){std::cout<<"q less than 1"<<std::endl;}
    else{
        q=q-1;
        double bb[16];
        int isignum=0;
        if(q>=2*nselau*nselau){
            TRANSPOSE_d4(b,bb);
            q=q-2*nselau*nselau;
            if(q>=nselau*nselau){
                q=q-nselau*nselau;
                isignum=1;
            }
        }
        else{
            for(int ii=0;ii<16;ii++){bb[ii]=b[ii];}
            if(q>=nselau*nselau){
                q=q-nselau*nselau;
                isignum=1;
            }
        }
        double btmp[16];
        int i,j;
        i=int(q/nselau);
        j=q % nselau;
        MATMAT_d4(&bsym[i*16],bb,btmp);
        MATMAT_d4(btmp,&bsym[j*16],bx);
        if(isignum==1){
            bx[3] = -bx[3];
            bx[7] = -bx[7];
            bx[11]= -bx[11];
            bx[12]= -bx[12];
            bx[13]= -bx[13];
            bx[14]= -bx[14];
        }
    }
}

void TJ::write_neighborInfo(std::string filename){
    std::ofstream outputFile(filename.c_str());
//    if(outputFile.is_open()){
//        for(unsigned int ii=0;ii<numGB;ii++){
//            for(int jj=0;jj<maxNeighbor;jj++){
//                outputFile << idxbs[ii*maxNeighbor+jj]<<"\t";
//            }
//            outputFile << std::endl;
//        }
//    }
    if(outputFile.is_open()){
        for(unsigned int ii=0;ii<numGB;ii++){
            for(int jj=0;jj<maxNeighbor*9;jj++){
                outputFile << Tranbs[ii*maxNeighbor*9+jj]<<"\t";
            }
            outputFile << std::endl;
        }
    }
    else{
        std::cout<<"Error opening: "<<filename<<std::endl;
    }
}

void TJ::write_idxcell(std::string filename){
    std::ofstream outputFile(filename.c_str());
    if(outputFile.is_open()){
        for(unsigned int ii=0;ii<numGB;ii++){
            for(int jj=0;jj<ksym;jj++){
                outputFile << idxcell[ii*ksym+jj]<<"\t";
            }
            outputFile << std::endl;
        }
    }
    else{
        std::cout<<"Error opening: "<<filename<<std::endl;
    }
}

void TJ::cal_idxcell(){
#pragma omp parallel for
   for(long int ii=0;ii<numGB;ii++){
       find_symeq_num(&bs[ii*16],ksym,
               &Qs[ii*ksym],&idxcell[ii*ksym]);
   } 
}

TJ::~TJ(){
//    delete idxcell;
//    delete Qs;
//    delete bs;
//    delete gcsym;
//    delete bsym;
}

void TJ::find_symeq_num(const double* b, int ksym,
        int* q, unsigned int* idxs ){
    bool testor;
    double bt[16];
    double bx[16];
    double btmp[16];
    double af[5];
    int isignum;
    TRANSPOSE_d4(b,bt);
    int nsymeq=-1;
    for(int ii=0;ii<nselau;ii++){
        for(int jj=0;jj<nselau;jj++){
            MATMAT_d4(&bsym[ii*16],b,btmp);
            MATMAT_d4(btmp,&bsym[jj*16],bx);
            matrix_to_param(bx,af,isignum);
            testor=subdomain.in_subdomain(af);
            if(testor){
                af[1]=1-cos(af[1]);
                af[3]=1-cos(af[3]);
                nsymeq+=1;
                if(nsymeq==ksym){return;}
                idxs[nsymeq]=subdomain.neuler_to_cell(af);
                q[nsymeq]=ii*nselau+jj+nselau*nselau*isignum+1; //q>=1 for any valid case
            }
            MATMAT_d4(&bsym[ii*16],bt,btmp);
            MATMAT_d4(btmp,&bsym[jj*16],bx);
            matrix_to_param(bx,af,isignum);
            testor=subdomain.in_subdomain(af);
            if(testor){
                af[1]=1-cos(af[1]);
                af[3]=1-cos(af[3]);
                nsymeq+=1;
                if(nsymeq==ksym){return;}
                idxs[nsymeq]=subdomain.neuler_to_cell(af);
                q[nsymeq]=ii*nselau+jj+2*nselau*nselau+nselau*nselau*isignum+1;
            }
        }
    }
}

void TJ::matrix_to_param(const double* b,
        double* fp, int &isignum){
    // first, calculate misorientation Euler
    double x = b[10];
    double sf = sqrt(1.0-std::min(x*x,1.0));
    if(sf>EPS){
        double tmp=-b[9]/sf;
        fp[0]=acoss(tmp);
        if(b[8]<0){fp[0]=2*PI-fp[0];}
        fp[1]=acoss(x);
        tmp=b[6]/sf;
        fp[2]=acoss(tmp);
        if(b[2]<0){fp[2]=2*PI-fp[2];}
    }
    else{
        fp[0]=b[0];
        fp[0]=acoss(fp[0]);
        if(b[1]<0){fp[0]=2*PI-fp[0];}
        fp[1]=acoss(x);
        fp[2]=0.0;
    }
    // then, calculate the plane normal
    fp[3]=acoss(b[11]);
    sf=sqrt(1.0-std::min(b[11]*b[11],1.0));
    if(sf>EPS){
        fp[4]=acoss(b[3]/sf);
        if(b[7]<0){fp[4]=2*PI-fp[4];}
    }
    else{
        fp[4]=acoss(b[3]);
        if(b[7]<0){fp[4]=2*PI-fp[4];} 
    }
    isignum=0;
    // ! change the direction to positive 'z'
    if(fp[3]>PI/2.0){
        fp[3]=PI-fp[3];
        fp[4]=fp[4]+PI;
        if(fp[4]> 2*PI){fp[4]=fp[4]-2*PI;}
        isignum=1;
    }
}

void TJ::fnormv(double* l){
    double norm=sqrt(l[0]*l[0]+l[1]*l[1]+l[2]*l[2]);
    for(int ii=0;ii<3;ii++){l[ii]=l[ii]/norm;}
}

void TJ::readSYM(std::string filename){
    std::ifstream inputFile(filename.c_str());
    if(inputFile.is_open()){
        std::string header;
        std::getline(inputFile,header);
        inputFile >> is;
        inputFile >> nse;
        inputFile >> nselau;
        gcsym = new double[nselau*3*3];
        bsym = new double[nselau*4*4](); // initialize to zero
        double tmp;
        for(int ii=0;ii<nselau;ii++){
            for(int jj=0;jj<3;jj++){
                for(int kk=0;kk<3;kk++){
                    inputFile >> tmp;
                    gcsym[ii*9+jj*3+kk]=tmp;
                    bsym[ii*16+jj*4+kk]=tmp;
                }
            }
            bsym[ii*16+15]=1;
        }
        inputFile.close();
    }
    else{
        std::cout<<"Error opening: "<<filename<<std::endl;
    }
}

void TJ::read_Config(std::string filename){
    std::ifstream inputFile(filename.c_str());
    double fmax[5];
    double tmp;
    int n[5];
    std::string TJfilename;
    std::string SYMfilename;
    if(inputFile.is_open()){
        std::string dump;
        std::getline(inputFile,dump);
        std::getline(inputFile,dump);
        std::getline(inputFile,dump);
        inputFile >> numTJ;
        std::getline(inputFile,dump);
        std::getline(inputFile,dump);
        std::getline(inputFile,dump);
        inputFile >> TJfilename;
        std::getline(inputFile,dump);
        std::getline(inputFile,dump);
        std::getline(inputFile,dump);
        inputFile >> maxNeighbor;
        std::getline(inputFile,dump);
        std::getline(inputFile,dump);
        std::getline(inputFile,dump);
        inputFile >> threshold;
        std::getline(inputFile,dump);
        std::getline(inputFile,dump);
        std::getline(inputFile,dump);
        std::getline(inputFile,dump);
        std::getline(inputFile,dump);
        inputFile >> SYMfilename;
        std::getline(inputFile,dump);
        std::getline(inputFile,dump);
        std::getline(inputFile,dump);
        for(int ii=0;ii<5;ii++){
            inputFile >> tmp;
            fmax[ii]=PI/tmp;
        }
        std::getline(inputFile,dump);
        std::getline(inputFile,dump);
        std::getline(inputFile,dump);
        for(int ii=0;ii<5;ii++){
            inputFile >> n[ii];
        }
        std::getline(inputFile,dump);
        std::getline(inputFile,dump);
        std::getline(inputFile,dump);
        inputFile >> ksym;
    }
    else{
        std::cout<<"Error opening: "<<filename<<std::endl;
    }
    numGB = numTJ*3;
    readTJ(TJfilename);
    readSYM(SYMfilename);
    subdomain.initialize(fmax,n,numGB,ksym);
    std::cout<<"TJfilename:"<<TJfilename<<std::endl;
    std::cout<<"SYMfilename:"<<SYMfilename<<std::endl;
    std::cout<<"numTJ:"<<numTJ<<std::endl;
    std::cout<<"ksym:"<<ksym<<std::endl;
}

TJ::TJ(std::string filename){
    read_Config(filename);
    idxcell= new unsigned int [numGB*ksym]();
    Qs = new int [numGB*ksym]();
}

void TJ::readTJ(std::string filename){
    std::ifstream inputFile(filename.c_str());
    if(inputFile.is_open()){
        bs= new double [numGB*4*4];
        double deg_rad=PI/180.0;
        double sl[3];
        double f[3];
        double sn1[3];
        double g1[9];
        double cn1[3];
        double gx1[9];
        double sn2[3];
        double g2[9];
        double cn2[3];
        double gx2[9];
        double sn3[3];
        double g3[9];
        double cn3[3];
        double gx3[9];
        unsigned int TJid;
        for(int ii=0;ii<numTJ;ii++){
            inputFile >> TJid;
            for(int jj=0;jj<3;jj++){inputFile>>sl[jj];}
            fnormv(sl);

                for(int jj=0;jj<3;jj++){
                    inputFile>>f[jj];
                    f[jj]*=deg_rad;
                }
                for(int jj=0;jj<3;jj++){inputFile>>sn1[jj];}
                fnormv(sn1);
                Euler2Mat(f,g1);
                for(int jj=0;jj<3;jj++){
                    inputFile>>f[jj];
                    f[jj]*=deg_rad;
                }
                for(int jj=0;jj<3;jj++){inputFile>>sn2[jj];}
                fnormv(sn2);
                Euler2Mat(f,g2);
                for(int jj=0;jj<3;jj++){
                    inputFile>>f[jj];
                    f[jj]*=deg_rad;
                }
                for(int jj=0;jj<3;jj++){inputFile>>sn3[jj];}
                fnormv(sn3);
                Euler2Mat(f,g3);

            MATVEC_d3(g2,sn1,cn1);
            MATVEC_d3(g3,sn2,cn2);
            MATVEC_d3(g1,sn3,cn3);
            
            double gtmp[9];
            TRANSPOSE_d3(g3,gtmp);
            MATMAT_d3(g2,gtmp,gx1);
            TRANSPOSE_d3(g1,gtmp);
            MATMAT_d3(g3,gtmp,gx2);
            TRANSPOSE_d3(g2,gtmp);
            MATMAT_d3(g1,gtmp,gx3);


                bconvert(gx1,cn1,&bs[ii*3*4*4+0*4*4]);
                bconvert(gx2,cn2,&bs[ii*3*4*4+1*4*4]);
                bconvert(gx3,cn3,&bs[ii*3*4*4+2*4*4]);
        }
        inputFile.close();
    }
    else{
        std::cout<<"Error opening: "<<filename<<std::endl;
    }
}


void TJ::find_neighbor(){
    std::cout<<"start finding neighbor..."<<std::endl;
    idxbs = new unsigned int [numGB*maxNeighbor]();
    disbs = new float [numGB*maxNeighbor]();
    Tranbs= new float [numGB*maxNeighbor*9]();
    NN = new int [numGB];

    subdomain.count_each_cell(idxcell,Qs);

#pragma omp parallel for
    for(long int kk=0;kk<numGB;kk++){
        unsigned int* cells=&idxcell[kk*ksym];
        int* theQs=&Qs[kk*ksym];
        double* b0=&bs[kk*16];
        double b1[16];
        double b2[16];
        int nn=0;
        for(int ii=0;ii<ksym;ii++){
            unsigned int cellid=cells[ii];
            int Q=theQs[ii];
            q2bx(Q,b0,b1);
            double M[9];
            normal_transformer(Q,b0,M);
            double Minv[9];
            TRANSPOSE_d3(M,Minv);
            for(int jj=0;jj<subdomain.num[cellid];jj++){
                bool firstMet=true;
                unsigned int bID=subdomain.ids[subdomain.cumnum[cellid]+jj];
                if(bID!=kk){
                    for(int ll=0;ll<nn;ll++){
                        firstMet=firstMet && (bID!=idxbs[kk*maxNeighbor+ll]);
                    }
                    if(firstMet){
                        int q=subdomain.qs[subdomain.cumnum[cellid]+jj];
                        q2bx(q,&bs[bID*16],b2);
                        double x= distance(b1,b2);
                        if(x<threshold){
                            idxbs[kk*maxNeighbor+nn]=bID;
                            disbs[kk*maxNeighbor+nn]=x;
                            double M2[9];
                            normal_transformer(q,&bs[bID*16],M2);
                            MATMAT_d3(Minv,M2,&Tranbs[kk*maxNeighbor*9+nn*9]);
                            nn+=1;
                        }
                    }
                }
                if(nn==maxNeighbor){break;}
            }
            if(nn==maxNeighbor){break;}
        }
        NN[kk]=nn;
    }
}

void TJ::Addrow(int row,int col,float val){
    rowA.push_back(row);
    colA.push_back(col);
    valA.push_back(val);
}

void TJ::writeA(){
    std::ofstream rowFile("rowA.binary",std::ios::binary);
    rowFile.write((char *)&rowA[0],rowA.size()*sizeof(int));
    rowFile.close();
    std::ofstream colFile("colA.binary",std::ios::binary);
    colFile.write((char *)&colA[0],colA.size()*sizeof(int));
    colFile.close();
    std::ofstream valFile("valA.binary",std::ios::binary);
    valFile.write((char *)&valA[0],valA.size()*sizeof(float));
    valFile.close();
}

void TJ::makeA(){
    int row_i=0;
    for(int kk=0;kk<numGB;kk++){
        for(int ii=0;ii<NN[kk];ii++){
            int idx=idxbs[kk*maxNeighbor+ii];
            float dis=disbs[kk*maxNeighbor+ii];
            float w=sqrt((1-dis/threshold+EPS)/NN[kk]);
            float* M=&Tranbs[kk*maxNeighbor*9+ii*9];
            Addrow(row_i,kk*3,-w);
            Addrow(row_i+1,kk*3+1,-w);
            Addrow(row_i+2,kk*3+2,-w);
            Addrow(row_i,idx*3,w*M[0]);
            Addrow(row_i,idx*3+1,w*M[1]);
            Addrow(row_i,idx*3+2,w*M[2]);
            Addrow(row_i+1,idx*3,w*M[3]);
            Addrow(row_i+1,idx*3+1,w*M[4]);
            Addrow(row_i+1,idx*3+2,w*M[5]);
            Addrow(row_i+2,idx*3,w*M[6]);
            Addrow(row_i+2,idx*3+1,w*M[7]);
            Addrow(row_i+2,idx*3+2,w*M[8]);
            row_i+=3;
        }
    }
}

void TJ::bconvert(const double* gx, const double* xn, double* b){
    for(int ii=0;ii<3;ii++){
        for(int jj=0;jj<3;jj++){
            b[ii*4+jj]=gx[ii*3+jj];
        }
    }
    for(int ii=0;ii<3;ii++){
        b[ii*4+3]=xn[ii];
    }
    double gtn[3];
    double gt[9];
    TRANSPOSE_d3(gx,gt);
    MATVEC_d3(gt,xn,gtn);
    for(int jj=0;jj<3;jj++){
        b[3*4+jj]=-gtn[jj];
    }
    b[15]=0.0;
}

