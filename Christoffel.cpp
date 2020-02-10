#include "headers/Calculus.h"
#include "headers/Metric.h"
#include "headers/Christoffel.h"

void InitC(double C[4][4][4]){
  for(int l=0;l<4;l++){
    for(int i=0;i<4;i++){
      for(int j=0;j<4;j++){
	C[l][i][j] = 0; 
      }
    }
  }
}

double SqrBrakets(MetricTensor g[4][4], double x[4], int i, int j, int k)
{
  double DiGkj = DPartial(i,g[k][j],x);
  double DjGki = DPartial(j,g[k][i],x);
  double DkGij = DPartial(k,g[i][j],x);

  return 0.5*(DiGkj + DjGki - DkGij);
}

double Gamma(MetricTensor g_ij[4][4], MetricTensor gij_[4][4], double x[4], int l, int i, int j)
{
  double ckGamma = 0;
  for(int k=0;k<4;k++){
      ckGamma += gij_[l][k](x)*SqrBrakets(g_ij,x,i,j,k);
  }
  /*if(abs(ckGamma) < dx)
    return 0;
    else*/
  return ckGamma;
}

void ChristoffelSymbol(double C[4][4][4],MetricTensor g_ij[4][4], MetricTensor gij_[4][4], double x[4])
{
  for(int l=0;l<4;l++){
    for(int i=0;i<4;i++){
      for(int j=0;j<4;j++){
	C[l][i][j] = Gamma(g_ij,gij_,x,l,i,j); 
      }
    }
  }
}

