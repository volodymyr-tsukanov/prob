/*
prob  Copyright (C) 2024  volodymyr-tsukanov
  This Source Code Form is subject to the terms of the Mozilla Public
  License, v. 2.0. If a copy of the MPL was not distributed with this
  file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/
#ifndef MINIMATH_H
#define MINIMATH_H
#include <cmath>


namespace VT::math {
static double logFactorial(int n){
    double logFact = 0;
    for (int i=1; i<=n; ++i){
        logFact += std::log(i);
    }
    return logFact;
}

static double binomialDistribution(int n, int k, double p){
    if(k==n) std::pow(p,n);  //probability of all successes
    if(k==0) return std::pow(1-p,n);  //probability of no successes
    if(k>n) return -1;

    double res = logFactorial(n) - logFactorial(k) - logFactorial(n - k);
    res += k * std::log(p) + (n-k) * std::log(1-p);
    return std::exp(res);
}
static double binomialDistributionPoisson(int n, int k, double p){
    if(k==n) std::pow(p,n);  //probability of all successes
    if(k==0) return std::pow(1-p,n);  //probability of no successes
    if(k>n) return -1;

    double a = n*p;
    double res = std::pow(a,k)*std::exp(-a) / std::tgamma(k+1); //tgamma(n) == (n-1)!
    return res;
}
static double binomialDistributionNormalApprox(int n, int k, double p){
    if(k==n) std::pow(p,n);  //probability of all successes
    if(k==0) return std::pow(1-p,n);  //probability of no successes
    if(k>n) return -1;

    double t = k-n*p, pw = t*t;
    t = 2*n*p*(1-p);
    double res = (1.0 / std::sqrt(t*M_PI)) * std::exp(-pw / t);
    return res;
}


static bool isInf(double n){
    return n == (1.0/0.0) || n == (-1.0/0.0);
}
static bool isNaN(double n) {
    return n != n;
}
}

#endif // MINIMATH_H
