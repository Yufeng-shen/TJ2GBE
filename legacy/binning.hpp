#ifndef binning_hpp
#define binning_hpp

#include<string>
#include<vector>

#define PI 3.14159265358979323846
#define EPS 0.000001 

//inline void MATMAT_d3(const double* A, const double* B,
//        double* res);
inline void MATMAT_d4(const double* A, const double* B,
        double* res);
inline void MATVEC_d3(const double* A, const double* V,
        double* res);
inline void Euler2Mat(const double* f, double* g);
inline void TRANSPOSE_d3(const double* A, double* res);
inline void TRANSPOSE_d4(const double* A, double* res);
inline double acoss(double ang);
inline double distance(const double* A, const double* B);

template <typename T1, typename T2, typename T3>
void MATMAT_d3(const T1 * A, const T2 * B,
        T3 * res){
    for(int ii=0;ii<3;ii++){
        for(int jj=0;jj<3;jj++){
            res[3*ii+jj]=0;
            for(int kk=0;kk<3;kk++){
                res[3*ii+jj]+=A[3*ii+kk]*B[3*kk+jj];
            }
        }
    }
}

class Subdomain{
    private:
        double fmax [5]; // subdomain max
        double afrange [5]; // subdomain range
        int    ndiv [5]; // number of cells in each dimension
        int totaln;
        unsigned int numGB;
        int ksym;
    public:
        unsigned int * num; // number of GB in cells
        unsigned int * cumnum;
        unsigned int * ids; // GB IDs in cells
        int * qs;

        void initialize(const double* fmax_, const int* n_,
                unsigned int numGB, int numSym);
        ~Subdomain();
        unsigned int neuler_to_cell(const double* af);
        bool in_subdomain(const double* af);
        void count_each_cell(const unsigned int* idxcell,const int* Qs);
};


class TJ{
    private:
        unsigned int * idxcell; // cell index of each GB
        int * Qs; // idx in equivalent set
        double* bs; // all GBs as 4x4 matrices
        double* gcsym;
        double* bsym; 
        unsigned int *idxbs; //index of neighbors
        float* disbs; //distance of neighbors
        float* Tranbs; //transformation for normal vector
        int * NN; //number of neighbors
        int maxNeighbor; //maximum number of neighbors
        float threshold; //threshold of distance for considering as neighbor
        int is,nse,nselau;
        unsigned int numTJ;
        unsigned int numGB;
        int ksym;
        std::vector <int> rowA;
        std::vector <int> colA;
        std::vector <float> valA;
        Subdomain subdomain;
        void fnormv(double* l); // l is length 3 vector
        void bconvert(const double* g, const double* xn,
                double* b);
        void find_symeq_num(const double* b, int ksym,
                int* q, unsigned int* idxs );
        void matrix_to_param(const double* b, 
                double* fp, int &isignum);
        void readTJ(std::string filename);
        void readSYM(std::string filename);
        void q2bx(int q, const double* b, double* bx);
        void normal_transformer(int q, const double* b, double* M);
        void Addrow(int row,int col,float val);
        void read_Config(std::string filename);
    public:
        TJ(std::string filename);
        ~TJ();
        void cal_idxcell();
        void write_idxcell(std::string filename);
        void write_neighborInfo(std::string filename);
        void write_binary();
        void find_neighbor();
        void makeA();
        void writeA();
};

#endif
