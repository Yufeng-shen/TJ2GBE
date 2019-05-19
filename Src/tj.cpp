#include "tj.hpp"

void TJ::cal_idxcell(){
#pragma omp parallel for
   for(long int ii=0;ii<numGB;ii++){
       find_symeq_num(bs[ii],cfg.ksym,
               &Qs[ii*cfg.ksym],&idxcell[ii*cfg.ksym]);
   } 
}

Matrix4<double> TJ::bconvert(const Matrix3<double> gx, const Vector3<double> xn){
    Matrix4<double> b;
    for(int ii=0;ii<3;ii++){
        for(int jj=0;jj<3;jj++){
            b.Set(ii,jj)=gx.Get(ii,jj);
        }
    }
    for(int ii=0;ii<3;ii++){
        b.Set(ii,3)=xn.Get(ii);
    }
    double gtn[3];
    double gt[9];
    TRANSPOSE_d3(gx,gt);
    MATVEC_d3(gt,xn,gtn);
    Vector3<double> gtn= gx.Transpose().dot(xn);
    for(int jj=0;jj<3;jj++){
        b.Set(3,jj)=-gtn.Get(jj);
    }
    b.Set(3,3)=0.0;
    return b;
}

void TJ::matrix_to_param(const Matrix4<double> b,
        double fp[5], int &isignum){
    // first, calculate misorientation Euler
    double x = b.Get(2,2);
    double sf = sqrt(1.0-std::min(x*x,1.0));
    if(sf>EPS){
        double tmp=-b.Get(2,1)/sf;
        fp[0]=acoss(tmp);
        if(b.Get(2,0)<0){fp[0]=2*PI-fp[0];}
        fp[1]=acoss(x);
        tmp=b.Get(1,2)/sf;
        fp[2]=acoss(tmp);
        if(b.Get(0,2)<0){fp[2]=2*PI-fp[2];}
    }
    else{
        fp[0]=b.Get(0,0);
        fp[0]=acoss(fp[0]);
        if(b.Get(0,1)<0){fp[0]=2*PI-fp[0];}
        fp[1]=acoss(x);
        fp[2]=0.0;
    }
    // then, calculate the plane normal
    fp[3]=acoss(b.Get(2,3));
    sf=sqrt(1.0-std::min(b.Get(2,3)*b.Get(2,3),1.0));
    if(sf>EPS){
        fp[4]=acoss(b.Get(0,3)/sf);
        if(b.Get(1,3)<0){fp[4]=2*PI-fp[4];}
    }
    else{
        fp[4]=acoss(b.Get(0,3));
        if(b.Get(1,3)<0){fp[4]=2*PI-fp[4];} 
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

void TJ::find_symeq_num(const Matrix4<double> b, int ksym,
        int* q, unsigned int* idxs ){
    bool testor;
    Matrix4<double> bt=b.Transpose();
    double af[5];
    int isignum;
    int nsymeq=-1;
    for(int ii=0;ii<nselau;ii++){
        for(int jj=0;jj<nselau;jj++){
            Matrix4<double> bx = bsym[ii].dot(b).dot(bsym[jj]);
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
            bx = bsym[ii].dot(bt).dot(bsym[jj]);
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


void TJ::Addrow(int row,int col,float val){
    rowA.push_back(row);
    colA.push_back(col);
    valA.push_back(val);
}


TJ::TJ(string filename){
    if(!cfg.InputConfigParameters(filename))
        exit(1);
    cfg.PrintFile();
    numGB = cfg.numTJ*3;
    readTJ(cfg.tripleJunctionFileName);
    readSYM(cfg.symmetryFileName);
    subdomain.initialize(cfg.fmax,cfg.n,numGB,cfg.ksym);
    idxcell= new unsigned int [numGB*ksym]();
    Qs = new int [numGB*ksym]();
}

void TJ::readTJ(string filename){
    ifstream inputFile(filename.c_str());
    if(inputFile.is_open()){
        bs.reserve(numGB);
        double deg_rad=PI/180.0;
        Vector3<double> sl;
        Vector3<double> f;
        Vector3<double> sn1;
        Matrix3<double> g1;
        Vector3<double> cn1;
        Matrix3<double> gx1;
        Vector3<double> sn2;
        Matrix3<double> g2;
        Vector3<double> cn2;
        Matrix3<double> gx2;
        Vector3<double> sn3;
        Matrix3<double> g3;
        Vector3<double> cn3;
        Matrix3<double> gx3;
        unsigned int TJid;
        for(int ii=0;ii<numTJ;ii++){
            inputFile >> TJid;
            for(int jj=0;jj<3;jj++){inputFile>>sl.Set(jj);}
            sl.Normalize();

                for(int jj=0;jj<3;jj++){
                    inputFile>>f.Set(jj);
                    f.Set(jj)*=deg_rad;
                }
                for(int jj=0;jj<3;jj++){inputFile>>sn1.Set(jj);}
                sn1.Normalize();
                g1=Euler2Mat(f);
                for(int jj=0;jj<3;jj++){
                    inputFile>>f.Set(jj);
                    f.Set(jj)*=deg_rad;
                }
                for(int jj=0;jj<3;jj++){inputFile>>sn2.Set(jj);}
                sn2.Normalize();
                g2=Euler2Mat(f);
                for(int jj=0;jj<3;jj++){
                    inputFile>>f.Set(jj);
                    f.Set(jj)*=deg_rad;
                }
                for(int jj=0;jj<3;jj++){inputFile>>sn3.Set(jj);}
                sn3.Normalize();
                g3=Euler2Mat(f);

            cn1=g2.dot(sn1);
            cn2=g3.dot(sn2);
            cn3=g1.dot(sn3);

            gx1=g2.dot(g3.Transpose());
            gx2=g3.dot(g1.Transpose());
            gx3=g1.dot(g2.Transpose());

            bs.push_back(bconvert(gx1,cn1));
            bs.push_back(bconvert(gx2,cn2));
            bs.push_back(bconvert(gx3,cn3));

        }
        inputFile.close();
    }
    else{
        std::cout<<"Error opening: "<<filename<<std::endl;
    }
}

void TJ::readSYM(string filename){
    ifstream inputFile(filename.c_str());
    if(inputFile.is_open()){
        string header;
        getline(inputFile,header);
        inputFile >> is;
        inputFile >> nse;
        inputFile >> nselau;
        gcsym.reserve(nselau);
        bsym.reserve(nselau);
        double tmp;
        for(int ii=0;ii<nselau;ii++){
            Matrix3 <double> gtmp;
            Matrix4 <double> btmp;
            btmp.Set(0,3)=0;
            btmp.Set(1,3)=0;
            btmp.Set(2,3)=0;
            btmp.Set(3,0)=0;
            btmp.Set(3,1)=0;
            btmp.Set(3,2)=0;
            btmp.Set(3,3)=1;
            for(int jj=0;jj<3;jj++){
                for(int kk=0;kk<3;kk++){
                    inputFile >> tmp;
                    gtmp.Set(jj,kk)=tmp;
                    btmp.Set(jj,kk)=tmp;
                }
            }
            gcsym.push_back(gtmp);
            bsym.push_back(btmp);
        }
        inputFile.close();
    }
    else{
        std::cout<<"Error opening: "<<filename<<std::endl;
    }
}
