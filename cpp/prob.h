#ifndef PROB_H
#define PROB_H
#include <iostream>
#include <initializer_list>
#include "miniMath.h"


namespace VT::prob {
enum P_precision{
    PRECISION_STRICT = 1, PRECISION_OPTIMAL = 6, PRECISION_FASTEST = 11
};

struct P_params {
    int n, k;
    double p;
};


class Possibility {
public:
    const int ERROR_INIT = -101;
    const double IGNORE_DBL = 1.0;


private:
    unsigned int calculation_count, recurrence_count;
    int k, prms_count;
    P_precision precision;
    P_params* prms;

protected:
    bool _dbg = false;

    double P(int n, int k, double p){   //P(X = k)
        if(p>=1) return IGNORE_DBL;
        ++calculation_count;
        double res = 0;
        if(n*p < 10 && n>100*p && (precision==PRECISION_FASTEST||precision==PRECISION_OPTIMAL)){
            if(_dbg) std::cout << "P";
            res = VT::math::binomialDistributionPoisson(n,k,p);
        } else if(n*p > 5 && n*(1-p) > 5 && precision==PRECISION_FASTEST){ //approx
            if(_dbg) std::cout << "N";
            res = VT::math::binomialDistributionNormalApprox(n,k,p);
        } else{
            if(_dbg) std::cout << "r";
            res = VT::math::binomialDistribution(n,k,p);
        }
        if(res > 1) res = IGNORE_DBL;
        return res;
    }
    double P_sum(int n, int k, double p){   //P(X >= k)
        if(k>n) return IGNORE_DBL;

        double res = 0;
        while(k<n){
            res += P(n, k++, p);
        }
        return res;
    }
    double P_combine_proc(int index, double** p){   //P(pass) = SUMM(k_1,k_2,...,k_n){ P(X_1=k_1)*...*P(X_n=k_n) }, where k_1+k_2+...+k_n >= k
        ++recurrence_count;
        double res = 0;
        if(index==0){  //tail => product calculation
            for(int k0=0; k0<=prms[0].n; k0++){
                if(prms[0].k+k0 >= k){
                    double product = p[0][k0];
                    if(_dbg) std::cout << k0<<":";
                    for(int c=1; c<prms_count; c++){
                        double a = p[c][prms[c].k];
                        product *= a;
                        if(_dbg) std::cout <<prms[c].k<<":";
                    }
                    if(_dbg) std::cout <<" "<<product << std::endl;
                    if(product<0 || product>1) product = 0;
                    res += product;
                }
            }
            return res;
        } else {  //head => summ calculation
            if(index == prms_count-1) prms[0].k = 0;
            for(int kn=0; kn<=prms[index].n; kn++){
                prms[0].k += kn;
                prms[index].k = kn;
                res += P_combine_proc(index-1, p);
                prms[0].k -= kn;
            }
            return res;
        }
    }
    double P_combine(){
        int c=0, n_max = prms[0].n;
        double** p = new double*[prms_count];
        for(c=1;c<prms_count;c++){  //calculate max n
            int n = prms[c].n;
            if(n > n_max) n_max = n;
        }
        for(c=0;c<prms_count;c++){ //pre-calculate partial possibilities
            p[c] = new double[n_max+1];
            if(_dbg) std::cout <<c<<"|\t";
            for(int i=0;i<=prms[c].n;i++){
                double a = P(prms[c].n,i,prms[c].p);
                if(VT::math::isInf(a) || VT::math::isNaN(a) || a<0) a = 0;
                p[c][i] = a;
                if(_dbg) std::cout <<a<<" ";
            }
            if(_dbg) std::cout << std::endl;
        }

        double res = P_combine_proc(prms_count-1, p);

        for(c=0;c<prms_count;c++) delete[] p[c];
        delete[] p;
        return res;
    }

public:
    Possibility(int n, int k, double p){
        precision = PRECISION_STRICT;
        prms_count = 1;
        prms = new P_params[prms_count];
        this->k = k;
        prms[0].n = n; prms[0].k = k; prms[0].p = p;
    }
    Possibility(int n, int k, double p, P_precision _precision) : Possibility(n,k,p){
        precision = _precision;
    }
    Possibility(int k, std::initializer_list<double> _prms){
        if(_prms.size()<2 || _prms.size()%2 != 0){
            prms_count = ERROR_INIT;
        } else {
            precision = PRECISION_OPTIMAL;
            prms_count = _prms.size()/2;
            prms = new P_params[prms_count];
            this->k = k;
            auto it = _prms.begin();
            for(int i=0; i<prms_count; i++){
                prms[i].n = (int) *(it++);
                prms[i].k = 0;
                prms[i].p = *(it++);
            }
        }
    }
    Possibility(const Possibility& otr){
        precision = otr.precision;
        prms_count = otr.prms_count;
        prms = new P_params[prms_count];
        for(int i=0; i<prms_count; i++){
            prms[i].n = otr.prms[i].n;
            prms[i].p = otr.prms[i].p;
        }
    }

    virtual ~Possibility(){
        if(prms_count > 0) delete[] prms;
    }

    Possibility& operator=(const Possibility& otr){
        if(this == &otr) return *this;

        delete[] prms;
        k = otr.k;
        prms_count = otr.prms_count;
        prms = new P_params[prms_count];
        for(int i=0; i<prms_count; i++){
            prms[i].n = otr.prms[i].n;
            prms[i].p = otr.prms[i].p;
        }

        return *this;
    }


    unsigned int getIterationCount(){
        return calculation_count+recurrence_count;
    }

    P_precision getPrecision(){
        return precision;
    }
    void setPrecision(P_precision _precision){
        precision = _precision;
    }


    double calculate(){
        calculation_count = recurrence_count = 0;
        double res = 0;
        if(prms_count <= 0){

        } else if(prms_count == 1){
            res = P_sum(prms[0].n,k,prms[0].p);
        } else{
            res = P_combine();
        }
        if(res>=1) res = 0.99999;
        return res;
    }

    void show(){
        if(prms_count <= 0) std::cout << "! not initialized\n";
        else{
            for(int i=0; i<prms_count; i++){
                std::cout << (i+1)<<": n="<<prms[i].n<<", k="<<prms[i].k<<", p="<<prms[i].p << std::endl;
            }
        }
    }
};
}

#endif // PROB_H
