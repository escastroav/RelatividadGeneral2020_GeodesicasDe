#include<iostream>
#include<functional>
#include<cmath>

#include "Metric.h"

double f01(MetricTensor g_ij[4][4], MetricTensor gij_[4][4], double x0[4], double x1[4], double x2[4], double t);
double f02(MetricTensor g_ij[4][4], MetricTensor gij_[4][4], double x0[4], double x1[4], double x2[4], double t);
double f03(MetricTensor g_ij[4][4], MetricTensor gij_[4][4], double x0[4], double x1[4], double x2[4], double t);
double f04(MetricTensor g_ij[4][4], MetricTensor gij_[4][4], double x0[4], double x1[4], double x2[4], double t);

double f11(MetricTensor g_ij[4][4], MetricTensor gij_[4][4], double x0[4], double x1[4], double x2[4], double t);
double f12(MetricTensor g_ij[4][4], MetricTensor gij_[4][4], double x0[4], double x1[4], double x2[4], double t);
double f13(MetricTensor g_ij[4][4], MetricTensor gij_[4][4], double x0[4], double x1[4], double x2[4], double t);
double f14(MetricTensor g_ij[4][4], MetricTensor gij_[4][4], double x0[4], double x1[4], double x2[4], double t);

void SolveGeodesics(MetricTensor g_ij[4][4], MetricTensor gij_[4][4], double (&x0)[4], double (&x1)[4], double (&x2)[4], double t);

