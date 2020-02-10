#include "headers/Geodesics.h"
#include "headers/Calculus.h"
#include "headers/Metric.h"
#include "headers/Christoffel.h"

double f01(MetricTensor g_ij[4][4], MetricTensor gij_[4][4], double x0[4], double x1[4], double x2[4], double t)
{
  return x2[0];
}
double f11(MetricTensor g_ij[4][4], MetricTensor gij_[4][4], double x0[4], double x1[4], double x2[4], double t)
{
  return x2[1];
}
double f21(MetricTensor g_ij[4][4], MetricTensor gij_[4][4], double x0[4], double x1[4], double x2[4], double t)
{
  return x2[2];
}
double f31(MetricTensor g_ij[4][4], MetricTensor gij_[4][4], double x0[4], double x1[4], double x2[4], double t)
{
  return x2[3];
}

double f02(MetricTensor g_ij[4][4], MetricTensor gij_[4][4], double x0[4], double x1[4], double x2[4], double t)
{
  double sum = 0;
  
  for(int i=0;i<4;i++){
    for(int j=0;j<4;j++){
      sum += Gamma(g_ij, gij_, x0, 0, i, j)*x1[i]*x1[j];
    }
  }
  return (-1)*sum;
}

double f12(MetricTensor g_ij[4][4], MetricTensor gij_[4][4], double x0[4], double x1[4], double x2[4], double t)
{
  double sum = 0;
  
  for(int i=0;i<4;i++){
    for(int j=0;j<4;j++){
      sum += Gamma(g_ij, gij_, x0, 1, i, j)*x1[i]*x1[j];
    }
  }
  return (-1)*sum;
}

double f22(MetricTensor g_ij[4][4], MetricTensor gij_[4][4], double x0[4], double x1[4], double x2[4], double t)
{
  double sum = 0;
  
  for(int i=0;i<4;i++){
    for(int j=0;j<4;j++){
      sum += Gamma(g_ij, gij_, x0, 2, i, j)*x1[i]*x1[j];
    }
  }
  return (-1)*sum;
}

double f32(MetricTensor g_ij[4][4], MetricTensor gij_[4][4], double x0[4], double x1[4], double x2[4], double t)
{
  double sum = 0;
  
  for(int i=0;i<4;i++){
    for(int j=0;j<4;j++){
      sum += Gamma(g_ij, gij_, x0, 3, i, j)*x1[i]*x1[j];
    }
  }
  return (-1)*sum;
}

void SolveGeodesics(MetricTensor g_ij[4][4], MetricTensor gij_[4][4],double (&x0)[4], double (&x1)[4], double (&x2)[4], double t)
{
  Equation f1[4] = {f01, f11, f21, f31};
  Equation f2[4] = {f02, f12, f22, f32};

  
cout << "tau" << "\t"
       << "r" << "\t"
       << "theta" << "\t"
       << "phi" << "\t"
       << "t" << "\t"
       << "dr" << "\t"
       << "dtheta" << "\t"
       << "dphi" << "\t"
       << "dt" << "\t"
       << "x" << "\t"
       << "y" << "\t"
       << endl;
  for(t=0; t<1e-1;t+=dt)
    {
      cout << t << "\t"
	   << x0[0] << "\t"
	   << x0[1] << "\t"
	   << x0[2] << "\t"
	   << x0[3] << "\t"
	   << x1[0] << "\t"
	   << x1[1] << "\t"
	   << x1[2] << "\t"
	   << x1[3] << "\t"
	   << x0[0]*cos(x0[2])*sin(x0[1]) << "\t"
	   << x0[0]*sin(x0[2])*sin(x0[1]) << "\t"
	   << endl;
      RK4(f1,f2,g_ij,gij_,x0,x1,x2,t);
    }
}

