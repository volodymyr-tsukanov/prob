/*
prob  Copyright (C) 2024  volodymyr-tsukanov
  This Source Code Form is subject to the terms of the Mozilla Public
  License, v. 2.0. If a copy of the MPL was not distributed with this
  file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/
#include <iostream>
#include <cmath>
#include "P.h"
#include "MiniMath.h"

using namespace std;
using namespace VT;


int main() {
    int n = 20;  //overall num
    int k = 11;  //num to pass
    double p4=0.25, p3=0.333, p2=0.5, p1=0.77997, p=0,t=0,res=0;  //partial probabilities
    prob::Possibility P (n, k, p4);


// 20: 1o4
    res = P.calculate();
    printf("P20 =\t%.5f=%.1fpercent  (%d iterations)\n",res,res*100,P.getIterationCount());
// 2: passed, 11: 1o3, 7: 1o2
    P = prob::Possibility(k-2, {11,p3, 7,p2});
    res = 2.0/n + P.calculate();
    printf("P2.11.7 =\t%.5f=%.1fpercent  (%d iterations)\n",res,res*100,P.getIterationCount());
// 2: passed, 11: 1o3, 7: 1o2
    P = prob::Possibility(k, {2,p1, 11,p3, 7,p2});
    res = P.calculate();
    printf("P2.11.7_t =\t%.5f=%.1fpercent  (%d iterations)\n",res,res*100,P.getIterationCount());
// 3: passed, 9: 1o2, 4: 1o3, 4: 1o4
    P = prob::Possibility(k, {3,p1, 9,p2, 4,p3, 4,p4});
    res = P.calculate();
    printf("P3.9.4.4 =\t%.5f=%.1fpercent  (%d iterations)\n",res,res*100,P.getIterationCount());
// 2: passed, 7: 1o2, 9: 1o3, 2: 1o4
    P = prob::Possibility(k, {2,p1, 7,p2, 9,p3, 2,p4});
    res = P.calculate();
    printf("P2.7.9.2 =\t%.5f=%.1fpercent  (%d iterations)\n",res,res*100,P.getIterationCount());
// 2: passed, 12: 1o2, 4: 1o3, 2: 1o4
    P = prob::Possibility(k, {2,p1, 12,p2, 4,p3, 2,p4});
    res = P.calculate();
    printf("P2.12.4.2 =\t%.5f=%.1fpercent  (%d iterations)\n",res,res*100,P.getIterationCount());


    // TESTS
    cout << "\n\n\n";
    n = 5; k = 5; p = 0.5;
    cout << "Test Case 1: (n = " << n << ", k = " << k << ", p = " << p << ", Expected: 0.03125)" << endl;
    res = prob::Possibility(n, k, p).calculate();  // Calculate res at the start of each test case
    t = math::binomialDistribution(n, k, p);
    cout << "\tr: " << t << "\t-> " << (t - res) << endl;
    t = math::binomialDistributionPoisson(n, k, p);
    cout << "\tP: " << t << "\t-> " << (t - res) << endl;
    t = math::binomialDistributionNormalApprox(n, k, p);
    cout << "\tn: " << t << "\t-> " << (t - res) << endl;

    n = 5; k = 0; p = 0.5;
    cout << "Test Case 2: (n = " << n << ", k = " << k << ", p = " << p << ", Expected: 0.03125)" << endl;
    res = prob::Possibility(n, k, p).calculate();  // Calculate res for this case
    t = math::binomialDistribution(n, k, p);
    cout << "\tr: " << t << "\t-> " << (t - res) << endl;
    t = math::binomialDistributionPoisson(n, k, p);
    cout << "\tp: " << t << "\t-> " << (t - res) << endl;
    t = math::binomialDistributionNormalApprox(n, k, p);
    cout << "\tn: " << t << "\t-> " << (t - res) << endl;

    n = 5; k = 2; p = 0.5;
    cout << "Test Case 3: (n = " << n << ", k = " << k << ", p = " << p << ", Expected: 0.3125)" << endl;
    res = prob::Possibility(n, k, p).calculate();  // Calculate res for this case
    t = math::binomialDistribution(n, k, p);
    cout << "\tr: " << t << "\t-> " << (t - res) << endl;
    t = math::binomialDistributionPoisson(n, k, p);
    cout << "\tp: " << t << "\t-> " << (t - res) << endl;
    t = math::binomialDistributionNormalApprox(n, k, p);
    cout << "\tn: " << t << "\t-> " << (t - res) << endl;

    n = 5; k = 6; p = 0.5;
    cout << "Test Case 4: (n = " << n << ", k = " << k << ", p = " << p << ", Expected: -1)" << endl;
    res = prob::Possibility(n, k, p).calculate();  // Calculate res for this case
    t = math::binomialDistribution(n, k, p);
    cout << "\tr: " << t << "\t-> " << (t - res) << endl;
    t = math::binomialDistributionPoisson(n, k, p);
    cout << "\tp: " << t << "\t-> " << (t - res) << endl;
    t = math::binomialDistributionNormalApprox(n, k, p);
    cout << "\tn: " << t << "\t-> " << (t - res) << endl;

    n = 10; k = 3; p = 0.7;
    cout << "Test Case 5: (n = " << n << ", k = " << k << ", p = " << p << ", Expected: 0.215233)" << endl;
    res = prob::Possibility(n, k, p).calculate();  // Calculate res for this case
    t = math::binomialDistribution(n, k, p);
    cout << "\tr: " << t << "\t-> " << (t - res) << endl;
    t = math::binomialDistributionPoisson(n, k, p);
    cout << "\tp: " << t << "\t-> " << (t - res) << endl;
    t = math::binomialDistributionNormalApprox(n, k, p);
    cout << "\tn: " << t << "\t-> " << (t - res) << endl;

    n = 0; k = 0; p = 0.7;
    cout << "Test Case 6: (n = " << n << ", k = " << k << ", p = " << p << ", Expected: 1)" << endl;
    res = prob::Possibility(n, k, p).calculate();  // Calculate res for this case
    t = math::binomialDistribution(n, k, p);
    cout << "\tr: " << t << "\t-> " << (t - res) << endl;
    t = math::binomialDistributionPoisson(n, k, p);
    cout << "\tp: " << t << "\t-> " << (t - res) << endl;
    t = math::binomialDistributionNormalApprox(n, k, p);
    cout << "\tn: " << t << "\t-> " << (t - res) << endl;

    n = 10; k = 0; p = 0.3;
    cout << "Test Case 7: (n = " << n << ", k = " << k << ", p = " << p << ", Expected: 0.00441)" << endl;
    res = prob::Possibility(n, k, p).calculate();  // Calculate res for this case
    t = math::binomialDistribution(n, k, p);
    cout << "\tr: " << t << "\t-> " << (t - res) << endl;
    t = math::binomialDistributionPoisson(n, k, p);
    cout << "\tp: " << t << "\t-> " << (t - res) << endl;
    t = math::binomialDistributionNormalApprox(n, k, p);
    cout << "\tn: " << t << "\t-> " << (t - res) << endl;

    n = 10; k = 8; p = 0.9;
    cout << "Test Case 8: (n = " << n << ", k = " << k << ", p = " << p << ", Expected: 0.193710244)" << endl;
    res = prob::Possibility(n, k, p).calculate();  // Calculate res for this case
    t = math::binomialDistribution(n, k, p);
    cout << "\tr: " << t << "\t-> " << (t - res) << endl;
    t = math::binomialDistributionPoisson(n, k, p);
    cout << "\tp: " << t << "\t-> " << (t - res) << endl;
    t = math::binomialDistributionNormalApprox(n, k, p);
    cout << "\tn: " << t << "\t-> " << (t - res) << endl;

    n = 10; k = 1; p = 0.1;
    cout << "Test Case 9: (n = " << n << ", k = " << k << ", p = " << p << ", Expected: 0.387420489)" << endl;
    res = prob::Possibility(n, k, p).calculate();  // Calculate res for this case
    t = math::binomialDistribution(n, k, p);
    cout << "\tr: " << t << "\t-> " << (t - res) << endl;
    t = math::binomialDistributionPoisson(n, k, p);
    cout << "\tp: " << t << "\t-> " << (t - res) << endl;
    t = math::binomialDistributionNormalApprox(n, k, p);
    cout << "\tn: " << t << "\t-> " << (t - res) << endl;

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

// GEN4
P20 =	0.00394=0.4percent  (9 iterations)
P2.11.7 =	0.53760=53.8percent  (63 iterations)
P2.11.7_t =	0.19068=19.1percent  (97 iterations)
P3.9.4.4 =	0.57638=57.6percent  (71 iterations)
P2.7.9.2 =	0.15689=15.7percent  (79 iterations)
P2.12.4.2 =	0.33801=33.8percent  (93 iterations)
*/
