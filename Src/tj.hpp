#ifndef TJ2GBE_TJ_HPP
#define TJ2GBE_TJ_HPP

#include "math.hpp"
#include "config.hpp"
#include "subdomain.hpp"
#include <vector>

using std::vector;
using std::string;
using std::ifstream;
using std::ofstream;

class TJ{
    private:
        Subdomain subdomain;
        Config cfg;
        unsigned int * idxcell; // cell index of each GB, shape: (numGB,ksym)
        int * Qs; // idx in equivalent set, shape: (numGB, ksym)
        vector< Matrix4<double> > bs;
        vector< Matrix3<double> > gcsym;
        vector< Matrix4<double> > bsym; 
        unsigned int *idxbs; //index of neighbors
        float* disbs; //distance of neighbors
        float* Tranbs; //transformation for normal vector
        vector<int> NN; //number of neighbors
        int is,nse,nselau;
        unsigned int numGB;
        vector <int> rowA;
        vector <int> colA;
        vector <float> valA;
        Matrix4<double> bconvert(const Matrix3<double> g, const Vector3<double> xn);
        void find_symeq_num(const Matrix4<double> b, int ksym,
                int* q, unsigned int* idxs );
        void matrix_to_param(const Matrix4<double> b, 
                double fp[5], int &isignum);
        void readTJ(string filename);
        void readSYM(string filename);
        Matrix4<double> q2bx(int q, const Matrix4<double> b);
        Matrix3<double> normal_transformer(int q, const Matrix4<double> b);
        void Addrow(int row,int col,float val);
    public:
        TJ(string filename);
        ~TJ();
        void cal_idxcell();
        void write_idxcell();
        void write_neighborInfo();
        void find_neighbor();
        void make_A();
        void write_A();
};

#endif
