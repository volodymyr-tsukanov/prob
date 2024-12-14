/*
prob  Copyright (C) 2024  volodymyr-tsukanov
  This Source Code Form is subject to the terms of the Mozilla Public
  License, v. 2.0. If a copy of the MPL was not distributed with this
  file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/
#include <iostream>
#include <cmath>
using namespace std;


struct P_params {
    int n, k;
    double p;
};


double logFactorial(int n){
    double logFact = 0;
    for (int i=1; i<=n; ++i){
        logFact += log(i);
    }
    return logFact;
}
bool isInf(double n){
    return n == (1.0/0.0) || n == (-1.0/0.0);
}
bool isNaN(double n) {
    return n != n;
}

double P(int n, int k, double p){
    if(k==n) pow(p,n);  //probability of all successes
    if(k==0) return pow(1-p,n);  //probability of no successes
    if(k>n) return -1;
    double res = logFactorial(n) - logFactorial(k) - logFactorial(n - k);
    res += k * log(p) + (n-k) * log(1-p);
    return exp(res);
}
double P_combine(int k, int n1, double p1, int n2, double p2){
    int k1=0, k2=0;
    double res = 0;
    for(int k1=0; k1<n1; k1++){
        for(int k2=0; k2<n2; k2++){
            if(k1+k2 != k) continue;
            double a = P(n1,k1,p1)*P(n2,k2,p2);
            if(isInf(a) || isNaN(a)) a = 0;
            res += a; //cout<<a<<" ";
        }
    }
    return res;
}
double P_combine_proc(const int& k, int index, double** p, P_params* prms, const int& _COUNT){
    if(prms[0].k > k) return 0;
    if(index==0){  //tail => product calculation
        for(int k0=0; k0<prms[0].n; k0++){
            if(prms[0].k+k0 == k){
                double product = p[0][k0];
                //cout <<k0<<":";
                for(int c=1; c<_COUNT; c++){
                    product *= p[c][prms[c].k];
                    //cout <<prms[c].k<<":";
                } //cout <<" "<<product << endl;
                return product;
            } //else cout <<"~"<<k0<<endl;
        }
        return 0;
    } else {  //head => summ calculation
        double res = 0;
        if(index == _COUNT-1) prms[0].k = 0;
        for(int kn=0; kn<prms[index].n; kn++){
            prms[0].k += kn;
            prms[index].k = kn;
            res += P_combine_proc(k,index-1, p,prms,_COUNT);
            prms[0].k -= kn;
        }
        return res;
    }
}
double P_combine(int k, P_params* prms, const int& _COUNT){
    if(_COUNT == 2) return P_combine(k, prms[0].n,prms[0].p, prms[1].n,prms[1].p);
    else if(_COUNT > 2){
        int i = 0, c = 0;
        double** p = new double*[_COUNT];
        for(c=0;c<_COUNT;c++){ //calculate partial possibilities
            p[c] = new double[k-1];
            for(i=1;i<k;i++){
                double a = P(prms[c].n,i,prms[c].p);
                if(isInf(a) || isNaN(a)) a = 0;
                p[c][i-1] = a; //cout <<a<<" ";
            } //cout<<endl;
        }

        //cout << "K=="<<k<<endl;
        double res = P_combine_proc(k,_COUNT-1, p,prms,_COUNT);

        for(int c=0;c<_COUNT;c++) delete p[c];
        delete[] p;
        return res;
    }
    return 0;
}


int main() {
    int n = 20;  //overall num
    int k = 11;  //num to pass
    double p4=0.25, p3=0.333, p2=0.5;  //partial probabilities
    P_params* pp = new P_params[5];

// 2: passed, 11: 1o3, 7: 1o2
    double res = 2.0/n + P_combine(k-2, 11,p3, 7,p2);
    printf("P2.11.7 = %.5f=%.1fpercent\n",res,res*100);
// 2: passed, 11: 1o3, 7: 1o2
    pp[0].n = 11; pp[0].p = p3;
    pp[1].n = 7; pp[1].p = p2;
    res = 2.0/n + P_combine(k-2, pp, 2);
    printf("P2.11.7_t = %.5f=%.1fpercent\n",res,res*100);


// 3: passed, 9: 1o2, 4: 1o3, 4: 1o4
    pp[0].n = 9; pp[0].p = p2;
    pp[1].n = 4; pp[1].p = p3;
    pp[2].n = 4; pp[2].p = p4;
    res = 3.0/n + P_combine(k-3, pp, 3);
    printf("P3.9.4.4 = %.5f=%.1fpercent\n",res,res*100);

// 2: passed, 7: 1o2, 9: 1o3, 2: 1o4
    pp[0].n = 7; pp[0].p = p2;
    pp[1].n = 9; pp[1].p = p3;
    pp[2].n = 2; pp[2].p = p4;
    res = 2.0/n + P_combine(k-2, pp, 3);
    printf("P2.7.9.2 = %.5f=%.1fpercent\n",res,res*100);

// 2: passed, 12: 1o2, 4: 1o3, 2: 1o4
    pp[0].n = 12; pp[0].p = p2;
    pp[1].n = 4; pp[1].p = p3;
    pp[2].n = 2; pp[2].p = p4;
    res = 2.0/n + P_combine(k-2, pp, 3);
    printf("P2.12.4.2 = %.5f=%.1fpercent\n",res,res*100);

    std::cout << "Test Case 1: " << P(5, 5, 0.5) << " (Expected: 0.03125)" << std::endl;
    // Test case 2
    std::cout << "Test Case 2: " << P(5, 0, 0.5) << " (Expected: 0.03125)" << std::endl;
    // Test case 3
    std::cout << "Test Case 3: " << P(5, 2, 0.5) << " (Expected: 0.3125)" << std::endl;
    // Test case 4
    std::cout << "Test Case 4: " << P(5, 6, 0.5) << " (Expected: -1)" << std::endl;
    // Test case 5
    std::cout << "Test Case 5: " << P(10, 3, 0.7) << " (Expected: ~0.215233)" << std::endl;

    delete pp;
    return 0;
}



/* OUTPUT
// GEN1
P2.11.7 = 0.35139=35.1percent
P3.9.4.4 = 0.24359=24.4percent
P2.7.9.2 = 0.14893=14.9percent
P2.12.4.2 = 0.20462=20.5percent

// GEN2
P2.11.7 = 0.47911=47.9percent
P3.9.4.4 = 0.63855=63.9percent
P2.7.9.2 = 0.25254=25.3percent
P2.12.4.2 = 0.83888=83.9percent
P2.11.7_t = 0.17170=17.2percent

// GEN3
P2.11.7 = 0.22571=22.6percent
P2.11.7_t = 0.22571=22.6percent
P3.9.4.4 = 0.23135=23.1percent
P2.7.9.2 = 0.10502=10.5percent
P2.12.4.2 = 0.10446=10.4percent
*/
