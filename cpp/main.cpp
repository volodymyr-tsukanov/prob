/*
prob  Copyright (C) 2024  volodymyr-tsukanov
  This Source Code Form is subject to the terms of the Mozilla Public
  License, v. 2.0. If a copy of the MPL was not distributed with this
  file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/
#include <iostream>
using namespace std;


long double factorial(int n){
    long double res = n;
    while(--n>1) res *= n;
    return res;
}
double power(double a, int p){
    while(--p>0) a += a;
    return a;
}

double P(int n, int k, double p){
    double res = 0;
    long double a = factorial(n) / factorial(k)*factorial(n-k);
    long double b = power(p,k);
    long double c = power((1.0-p),n-k);
    res = a*b*c;
    cout << a<<" "<<b<<" "<<c<<" ";
    return res;
}


int main() {
    int n = 20;  //overall num
    int k = 11;  //num to pass
    double p=0.25;  //probability of single pass

    cout << "P="<<P(n,k,p) << endl;
    return 0;
}
