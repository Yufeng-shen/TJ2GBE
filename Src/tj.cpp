#include "tj.hpp"

void TJ::cal_idxcell(){
#pragma omp parallel for
   for(long int ii=0;ii<numGB;ii++){
       find_symeq_num(bs[ii],cfg.ksym,
               &Qs[ii*cfg.ksym],&idxcell[ii*cfg.ksym]);
   } 
}

void TJ::find_neighbor(){
    std::cout<<"start finding neighbor..."<<std::endl;
    idxbs = new unsigned int [numGB*cfg.maxNeighbor]();
    disbs = new float [numGB*cfg.maxNeighbor]();
    Tranbs= new float [numGB*cfg.maxNeighbor*9]();
    NN.reserve(numGB);

    subdomain.count_each_cell(idxcell,Qs);

#pragma omp parallel for
    for(long int kk=0;kk<numGB;kk++){
        unsigned int* cells=&idxcell[kk*cfg.ksym];
        int* theQs=&Qs[kk*cfg.ksym];
        Matrix4<double> b0=bs[kk];
        int nn=0;
        for(int ii=0;ii<cfg.ksym;ii++){
            unsigned int cellid=cells[ii];
            int Q=theQs[ii];
            Matrix4<double> b1=q2bx(Q,b0);
            Matrix3<double> M=normal_transformer(Q,b0);
            Matrix3<double> Minv=M.Transpose();
            for(int jj=0;jj<subdomain.num[cellid];jj++){
                bool firstMet=true;
                unsigned int bID=subdomain.ids[subdomain.cumnum[cellid]+jj];
                if(bID!=kk){
                    for(int ll=0;ll<nn;ll++){
                        firstMet=firstMet && (bID!=idxbs[kk*cfg.maxNeighbor+ll]);
                    }
                    if(firstMet){
                        int q=subdomain.qs[subdomain.cumnum[cellid]+jj];
                        Matrix4<double> b2=q2bx(q,bs[bID]);
                        double x= distance(b1,b2);
                        if(x<threshold){
                            idxbs[kk*cfg.maxNeighbor+nn]=bID;
                            disbs[kk*cfg.maxNeighbor+nn]=x;
                            Matrix3<double> M2=normal_transformer(q,bs[bID]);
                            Matrix3<float> M=Minv.dot(M2);
                            for(int iii=0;iii<3;iii++)
                                for(int jjj=0;jjj<3;jjj++)
                                    Tranbs[kk*cfg.maxNeighbor*9+nn*9+iii*3+jjj]=M(iii,jjj);
                            nn+=1;
                        }
                    }
                }
                if(nn==cfg.maxNeighbor){break;}
            }
            if(nn==cfg.maxNeighbor){break;}
        }
        NN[kk]=nn;
    }
}

void TJ::make_A(){
    int row_i=0;
    for(int kk=0;kk<numGB;kk++){
        for(int ii=0;ii<NN[kk];ii++){
            int idx=idxbs[kk*cfg.maxNeighbor+ii];
            float dis=disbs[kk*cfg.maxNeighbor+ii];
            float w=sqrt((1-dis/cfg.threshold+EPS)/NN[kk]);
            float* M=&Tranbs[kk*cfg.maxNeighbor*9+ii*9];
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

void TJ::write_idxcell(){
    ofstream outputFile("idxcell.txt");
    if(outputFile.is_open()){
        for(unsigned int ii=0;ii<numGB;ii++){
            for(int jj=0;jj<cfg.ksym;jj++){
                outputFile << idxcell[ii*cfg.ksym+jj]<<"\t";
            }
            outputFile << std::endl;
        }
    }
    else{
        std::cout<<"Writing idxcell.txt failed"<<std::endl;
    }
}

void TJ::write_neighborInfo(){
    ofstream outputFile("idxbs.txt");
    if(outputFile.is_open()){
        for(unsigned int ii=0;ii<numGB;ii++){
            for(int jj=0;jj<cfg.maxNeighbor;jj++){
                outputFile << idxbs[ii*cfg.maxNeighbor+jj]<<"\t";
            }
            outputFile << std::endl;
        }
    }
    else{
        std::cout<<"Writing idxbs.txt failed"<<std::endl;
    }
    ofstream outputFile2("disbs.txt");
    if(outputFile2.is_open()){
        for(unsigned int ii=0;ii<numGB;ii++){
            for(int jj=0;jj<cfg.maxNeighbor;jj++){
                outputFile2 << disbs[ii*cfg.maxNeighbor+jj]<<"\t";
            }
            outputFile2 << std::endl;
        }
    }
    else{
        std::cout<<"Writing disbs.txt failed"<<std::endl;
    }
    ofstream outputFile3("NN.txt");
    if(outputFile3.is_open()){
        for(unsigned int ii=0;ii<numGB;ii++){
                outputFile3 << NN[ii]<<std::endl;
        }
    }
    else{
        std::cout<<"Writing NN.txt failed"<<std::endl;
    }
}

void TJ::write_A(){
    ofstream rowFile("rowA.binary",std::ios::binary);
    rowFile.write((char *)&rowA[0],rowA.size()*sizeof(int));
    rowFile.close();
    ofstream colFile("colA.binary",std::ios::binary);
    colFile.write((char *)&colA[0],colA.size()*sizeof(int));
    colFile.close();
    ofstream valFile("valA.binary",std::ios::binary);
    valFile.write((char *)&valA[0],valA.size()*sizeof(float));
    valFile.close();
}


Matrix3<double> TJ::normal_transformer(int q, const Matrix4<double> b){
    Matrix3<double> M;
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
        M=gcsym[i];
    }
    else if((isignum)&&(!istranspose)){
        M=-gcsym[i]
    }
    else if((!isignum)&&(istranspose)){
        Matrix3<double> tmp;
        for(int iii=0;iii<3;iii++){
            for(int jjj=0;jjj<3;jjj++){
                tmp.Set(iii,jjj)=-b.Get(jjj,iii);
            }
        }
        M=gcsym[i].dot(tmp);
    }
    else{
        Matrix3<double> tmp;
        for(int iii=0;iii<3;iii++){
            for(int jjj=0;jjj<3;jjj++){
                tmp.Set(iii,jjj)=b.Get(jjj,iii);
            }
        }
        M=gcsym[i].dot(tmp);
    }
}

Matrix4<double> TJ::q2bx(int q, const Matrix4<double> b){
    Matrix4<double> bx;
    if(q<1){
        std::cout<<"Error: q less than 1"<<std::endl;
        return bx;
    }
    else{
        q=q-1;
        int isignum=0;
        if(q>=2*nselau*nselau){
            Matrix4<double> bb= b.Transpose();
            q=q-2*nselau*nselau;
            if(q>=nselau*nselau){
                q=q-nselau*nselau;
                isignum=1;
            }
        }
        else{
            Matrix4<double> bb=b;
            if(q>=nselau*nselau){
                q=q-nselau*nselau;
                isignum=1;
            }
        }
        int i,j;
        i=int(q/nselau);
        j=q % nselau;
        bx=bsym[i].dot(bb).dot(bsym[j]);
        if(isignum==1){
            bx.Set(0,3) = -bx.Get(0,3);
            bx.Set(1,3) = -bx.Get(1,3);
            bx.Set(2,3) = -bx.Get(2,3);
            bx.Set(3,0) = -bx.Get(3,0);
            bx.Set(3,1) = -bx.Get(3,1);
            bx.Set(3,2) = -bx.Get(3,2);
        }
        return bx;
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
    idxcell= new unsigned int [numGB*cfg.ksym]();
    Qs = new int [numGB*cfg.ksym]();
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

TJ::~TJ(){
    delete idxcell;
    delete Qs;
    delete idxbs;
    delete disbs;
    delete Tranbs;
}
