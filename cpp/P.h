/*
prob  Copyright (C) 2024  volodymyr-tsukanov
  This Source Code Form is subject to the terms of the Mozilla Public
  License, v. 2.0. If a copy of the MPL was not distributed with this
  file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/
#ifndef P_H
#define P_H
#include <iostream>
#include <initializer_list>
#include "MiniMath.h"


namespace VT::prob {
struct P_params {
    int n, k;
    double p;
};


class Possibility {
public:
    const int ERROR_INIT = -101;
    const double IGNORE_DBL = 1.0;


private:
    unsigned int iteration_count;
    int k, prms_count;
    P_params* prms;

protected:
    double P(int n, int k, double p){
        ++iteration_count;
        if(n*p < 10 && n>100*p){
            //std::cout << "P";
            return VT::math::binomialDistributionPoisson(n,k,p);
        } else if(n*p > 5 && n*(1-p) > 5){ //approx
            //std::cout << "n";
            return VT::math::binomialDistributionNormalApprox(n,k,p);
        } else{
            //std::cout << "r";
            return VT::math::binomialDistribution(n,k,p);
        }
    }
    double P_sum(int n, int k, double p){
        if(k>n) return IGNORE_DBL;

        double res = 0;
        while(k<n){
            res += P(n, k++, p);
        }
        return res;
    }
    double P_combine(int n1, double p1, int n2, double p2){
        int k1=0, k2=0;
        double res = 0;
        for(int k1=0; k1<n1; k1++){
            for(int k2=0; k2<n2; k2++){
                if(k1+k2 != k) continue;
                double a = P_sum(n1,k1,p1)*P_sum(n2,k2,p2);
                if(VT::math::isInf(a) || VT::math::isNaN(a)) a = 0;
                res += a; //cout<<a<<" ";
            }
        }
        return res;
    }
    double P_combine_proc(int index, double** p){
        if(prms[0].k > k) return 0;
        if(index==0){  //tail => product calculation
            for(int k0=0; k0<prms[0].n; k0++){
                if(prms[0].k+k0 == k){
                    double product = 1.0;
                    for(int c=0; c<prms_count; c++){
                        double a = p[c][prms[c].k];
                        product *= a;
                        //cout <<prms[c].k<<":";
                    } //cout <<" "<<product << endl;
                    return product;
                } //else cout <<"~"<<k0<<endl;
            }
            return 0;
        } else {  //head => summ calculation
            double res = 0;
            if(index == prms_count-1) prms[0].k = 0;
            for(int kn=0; kn<prms[index].n; kn++){
                prms[0].k += kn;
                prms[index].k = kn;
                res += P_combine_proc(index-1, p);
                prms[0].k -= kn;
            }
            return res;
        }
    }
    double P_combine(){
        int c = 0;
        double** p = new double*[prms_count];
        for(c=0;c<prms_count;c++){ //calculate partial possibilities
            p[c] = new double[k];
            for(int i=0;i<k;i++){
                double a = P_sum(prms[c].n,i,prms[c].p);
                if(VT::math::isInf(a) || VT::math::isNaN(a)) a = 0;
                p[c][i] = a; //cout <<a<<" ";
            } //cout<<endl;
        }

        //cout << "K=="<<k<<endl;
        double res = P_combine_proc(prms_count-1, p);

        for(c=0;c<prms_count;c++) delete[] p[c];
        delete[] p;
        return res;
    }

public:
    Possibility(int n, int k, double p){
        iteration_count = 0;
        prms_count = 1;
        prms = new P_params[prms_count];
        this->k = k;
        prms[0].n = n; prms[0],k = k; prms[0].p = p;
    }
    Possibility(int k, std::initializer_list<double> _prms){
        if(_prms.size()<2 || _prms.size()%2 != 0){
            prms_count = ERROR_INIT;
        } else {
            iteration_count = 0;
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
        return iteration_count;
    }


    double calculate(){
        iteration_count = 0;
        if(prms_count <= 0){

        } else if(prms_count == 1){
            return P_sum(prms[0].n,k,prms[0].p);
        } else if(prms_count == 2){
            return P_combine(prms[0].n,prms[0].p, prms[1].n,prms[1].p);
        } else {
            return P_combine();
        }
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

#endif // P_H
