#include "headers/Metric.h"
#include "headers/Constants.h"

double zero(double x[4]){return 0;};

void InitMetric(MetricTensor g[4][4]){
  for(int i=0;i<4;i++){
    for(int j=0;j<4;j++){
	g[i][j] = zero;
    }
  }
}

//Schwarzschild covariant
double g_00(double x[4]){return -1/(1-2*Rs/x[0]);}
double g_11(double x[4]){return (-1)*pow(x[0],2.0);}
double g_22(double x[4]){return (-1)*pow(x[0]*sin(x[1]),2.0);}
double g_33(double x[4]){return (1-2*Rs/x[0]);}

//Schwarzschild contravariant
double g00_(double x[4]){return (-1)*(1-2*Rs/x[0]);}
double g11_(double x[4]){return (-1)*pow(x[0],-2.0);}
double g22_(double x[4]){return (-1)*pow(x[0]*sin(x[1]),-2.0);}
double g33_(double x[4]){return 1/(1-2*Rs/x[0]);}
