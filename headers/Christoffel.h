#include<iostream>
#include<functional>
#include<cmath>


using namespace std;

void InitC(double C[4][4][4]);
double SqrBrakets(MetricTensor g[4][4], double x[4], int i, int j, int k);
double Gamma(MetricTensor g_ij[4][4], MetricTensor gij_[4][4], double x[4], int l, int i, int j);
void ChristoffelSymbol(double C[4][4][4],MetricTensor g_ij[4][4], MetricTensor gij_[4][4], double x[4]);




