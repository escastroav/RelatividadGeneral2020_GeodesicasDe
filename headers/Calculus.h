#include<iostream>
#include<functional>
#include<cmath>

#include "Metric.h"

typedef double (*Equation) (MetricTensor g_ij[4][4], MetricTensor gij_[4][4], double * x0, double * x1, double * x2, double t);

using namespace std;

const double dx=1.3e-8;
const double dt=1.3e-8;

double DPartial(int i, function<double (double *)> g,double x[4]);
void RK4(Equation f1[4],Equation f2[4], MetricTensor g_ij[4][4], MetricTensor gij_[4][4], double (&x0)[4], double (&x1)[4], double (&x2)[4], double t);

