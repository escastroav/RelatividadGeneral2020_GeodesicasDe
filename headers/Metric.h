#include<iostream>
#include<functional>
#include<cmath>

using namespace std;

typedef double (*MetricTensor) (double * x);
double zero(double x[4]);

void InitMetric(MetricTensor g[4][4]);

//Covariant
double g_00(double x[4]);
double g_11(double x[4]);
double g_22(double x[4]);
double g_33(double x[4]);

//Contravariant
double g00_(double x[4]);
double g11_(double x[4]);
double g22_(double x[4]);
double g33_(double x[4]);
