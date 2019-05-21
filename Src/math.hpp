#ifndef TJ2GBE_MATH_HPP
#define TJ2GBE_MATH_HPP

#include <cmath>
#include <iostream>

#ifndef PI
#define PI 3.14159265358979323846
#endif

#ifndef EPS
#define EPS 0.000001 
#endif

template <typename _Scalar>
class Vector3{
    private:
        _Scalar v[3];
    public:
        Vector3(){}
        template<typename _Scalar2>
        Vector3(const Vector3 < _Scalar2> &rhs){
            v[0]=rhs.Get(0);
            v[1]=rhs.Get(1);
            v[2]=rhs.Get(2);
        }
        inline _Scalar & Set(int i){return v[i];}
        inline _Scalar Get(int i) const {return v[i];}
        inline void Normalize(){
            _Scalar norm=sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
            for(int ii=0;ii<3;ii++)
                v[ii]=v[ii]/norm;
        }
};


template <typename _Scalar>
class Matrix3{
    private:
        _Scalar m[3][3];
    public:
        Matrix3(){}
        template<typename _Scalar2>
        Matrix3(const Matrix3<_Scalar2> &rhs){
            for(int ii=0;ii<3;ii++){
                for(int jj=0;jj<3;jj++){
                    m[ii][jj]=rhs.Get(ii,jj);
                }
            }
        }    
        template <typename _Scalar2>
        inline void operator=( Matrix3<_Scalar2> &rhs){
            for(int ii=0;ii<3;ii++){
                for(int jj=0;jj<3;jj++){
                    m[ii][jj]=rhs.Get(ii,jj);
                }
            }
        }
        inline Matrix3< _Scalar> operator-() const {
            Matrix3<_Scalar> res;
            for(int ii=0;ii<3;ii++){
                for(int jj=0;jj<3;jj++){
                    res.Set(ii,jj)=-m[ii][jj];
                }
            }
            return res;
        }

        inline _Scalar & Set(int i,int j){return m[i][j];}
        inline _Scalar Get(int i,int j) const {return m[i][j];}
        inline Matrix3 < _Scalar > Transpose() const{
            Matrix3 < _Scalar > res;
            for(int ii=0;ii<3;ii++){
                for(int jj=0;jj<3;jj++)
                    res.Set(ii,jj)=m[jj][ii];
            }
            return res;
        }
        template <typename _Scalar2>
        inline Vector3 < _Scalar> dot(const Vector3 <_Scalar2> & rhsV) const{
            Vector3 < _Scalar> res;
            for(int ii=0; ii<3;ii++){
                res.Set(ii)=0;
                for(int jj=0;jj<3;jj++){
                    res.Set(ii)+=m[ii][jj]*rhsV.Get(jj);
                }
            }
            return res;
        }
        template <typename _Scalar2>
        inline Matrix3 < _Scalar > dot(const Matrix3 < _Scalar2> & rhsM) const{
            Matrix3 < _Scalar > res;
            for(int ii=0;ii<3;ii++){
                for(int jj=0;jj<3;jj++){
                    res.Set(ii,jj)=0;
                    for(int kk=0;kk<3;kk++){
                        res.Set(ii,jj)+=m[ii][kk]*rhsM.Get(kk,jj);
                    }
                }
            }
            return res;
        }
        void Print() const{
            std::cout<<m[0][0]<<"\t"<<m[0][1]<<"\t"<<m[0][2]<<std::endl
            <<m[1][0]<<"\t"<<m[1][1]<<"\t"<<m[1][2]<<std::endl
            <<m[2][0]<<"\t"<<m[2][1]<<"\t"<<m[2][2]<<std::endl;
            return;
        }
};

template <typename _Scalar>
class Matrix4{
    private:
        _Scalar b[4][4];
    public:
        Matrix4(){}
        template <typename _Scalar2>
        Matrix4(const Matrix4<_Scalar2> &rhs){
            for(int ii=0;ii<4;ii++){
                for(int jj=0;jj<4;jj++){
                    b[ii][jj]=rhs.Get(ii,jj);
                }
            }
        }   
        template <typename _Scalar2>
        inline void operator=( Matrix4<_Scalar2> &rhs){
            for(int ii=0;ii<4;ii++){
                for(int jj=0;jj<4;jj++){
                    b[ii][jj]=rhs.Get(ii,jj);
                }
            }
        }
        inline _Scalar & Set(int i,int j){return b[i][j];}
        inline _Scalar Get(int i,int j) const {return b[i][j];}
        inline Matrix4 < _Scalar > Transpose() const {
            Matrix4 < _Scalar > res;
            for(int ii=0;ii<4;ii++){
                for(int jj=0;jj<4;jj++)
                    res.Set(ii,jj)=b[jj][ii];
            }
            return res;
        }
        template <typename _Scalar2>
        inline Matrix4 <_Scalar> dot(const Matrix4 < _Scalar2> & rhsB) const {
            Matrix4 < _Scalar > res;
            for(int ii=0;ii<4;ii++){
                for(int jj=0;jj<4;jj++){
                    res.Set(ii,jj)=0;
                    for(int kk=0;kk<4;kk++){
                        res.Set(ii,jj)+=b[ii][kk]*rhsB.Get(kk,jj);
                    }
                }
            }
            return res;
        }
        void Print() const{
            std::cout<<b[0][0]<<"\t"<<b[0][1]<<"\t"<<b[0][2]<<"\t"<<b[0][3]<<std::endl
            <<b[1][0]<<"\t"<<b[1][1]<<"\t"<<b[1][2]<<"\t"<<b[1][3]<<std::endl
            <<b[2][0]<<"\t"<<b[2][1]<<"\t"<<b[2][2]<<"\t"<<b[2][3]<<std::endl
            <<b[3][0]<<"\t"<<b[3][1]<<"\t"<<b[3][2]<<"\t"<<b[3][3]<<std::endl;
            return;
        }
};

inline Matrix3< double> Euler2Mat(const Vector3< double> & f){
    Matrix3< double > g;
    double sf=sin(f.Get(1));
    double cf=cos(f.Get(1));
    double sf1=sin(f.Get(0));
    double cf1=cos(f.Get(0));
    double sf2=sin(f.Get(2));
    double cf2=cos(f.Get(2));
    g.Set(0,0)=cf1*cf2-sf1*sf2*cf;
  
    g.Set(0,1)=sf1*cf2+cf1*sf2*cf;
  
    g.Set(0,2)=sf2*sf;
   
    g.Set(1,0)=-cf1*sf2-sf1*cf2*cf;
  
    g.Set(1,1)=-sf1*sf2+cf1*cf2*cf;
  
    g.Set(1,2)=cf2*sf;
   
    g.Set(2,0)=sf1*sf;
  
    g.Set(2,1)=-cf1*sf;
  
    g.Set(2,2)=cf;

    return g;
}

inline double acoss(double ang){
    if(ang>1) return 0.0;
    if(ang<-1) return PI;
    return acos(ang);
}

template <typename _Scalar>
inline _Scalar distance(const Matrix4< _Scalar > &A, const Matrix4< _Scalar > &B){
    _Scalar res=0;
    for(int ii=0;ii<4;ii++){
        for(int jj=0;jj<4;jj++)
            res += (A.Get(ii,jj)-B.Get(ii,jj))*(A.Get(ii,jj)-B.Get(ii,jj));
    }
    return res;
}

#endif
