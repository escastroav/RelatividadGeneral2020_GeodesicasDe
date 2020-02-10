#include "headers/Calculus.h"
#include "headers/Metric.h"
#include "headers/Christoffel.h"
#include "headers/Constants.h"
#include "headers/Geodesics.h"


double Potential(double x[4], double epsilon, double L)
{
  return (0.5*epsilon - 0.5*Rs*epsilon/x[0] + L*L*0.5/(x[0]*x[0]) - 0.5*Rs*L*L/(x[0]*x[0]*x[0]));
}

int main(void){
  cout.precision(4);
  double t = 0;

  double r0 = 5*Rs;
  double T0 = 5e-8;
  double L0 = 20*M_PI;
  double E0 = 1.0;

  MetricTensor g_ij[4][4]; InitMetric(g_ij);
  MetricTensor gij_[4][4]; InitMetric(gij_);
  
  double x0[4] = {r0, 0.5*M_PI, 0, 0};
  double x1[4] = {0.5*E0*E0 - Potential(x0, 1, L0), 0, L0/(r0*r0), E0 / g_00(x0)}; 
  double x2[4] = {-Rs*0.5/(r0*r0),0,0,0};  

  g_ij[0][0] = g_00;
  g_ij[1][1] = g_11;
  g_ij[2][2] = g_22;
  g_ij[3][3] = g_33;

  gij_[0][0] = g00_;
  gij_[1][1] = g11_;
  gij_[2][2] = g22_;
  gij_[3][3] = g33_;

  SolveGeodesics(g_ij, gij_, x0, x1, x2, t);

  
  
  return 0;
}
