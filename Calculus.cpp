#include "headers/Calculus.h"
#include "headers/Metric.h"
#include "headers/Geodesics.h"


double DPartial(int i, function<double (double *)> g,double x[4]){
  /* if(g(x) == 0)
     return 0;*/
  
  double d1,d2,d3,d4;
  
  x[i] -= 2*dx;
  d1 = g(x);
  x[i] += dx;
  d2 = g(x);
  x[i] += 2*dx;
  d3 = g(x);
  x[i] += dx;
  d4 = g(x);
  x[i] -= 2*dx;
  
  return (d1-8*d2+8*d3-d4)/(12*dx);
}

void RK4(Equation f1[4], Equation f2[4], MetricTensor g_ij[4][4], MetricTensor gij_[4][4], double (&x0)[4], double (&x1)[4], double (&x2)[4], double t)
{
  double dx01 [4],dx02 [4],dx03 [4],dx04 [4]; 
  double dx11 [4],dx12 [4],dx13 [4],dx14 [4]; 
  double dx21 [4],dx22 [4],dx23 [4],dx24 [4];
  
  for(int i=0;i<4;i++) dx01[i]=f1[i](g_ij,gij_,x0,x1,x2,t)*dt;
for(int i=0;i<4;i++)   dx11[i]=f1[i](g_ij,gij_,x0,x1,x2,t)*dt;
for(int i=0;i<4;i++)    dx21[i]=f2[i](g_ij,gij_,x0,x1,x2,t)*dt;
  
for(int i=0;i<4;i++)   x0[i] += 0.5*dx01[i]; 
for(int i=0;i<4;i++)   x1[i] += 0.5*dx11[i]; 
for(int i=0;i<4;i++)   x2[i] += 0.5*dx21[i]; 
  
for(int i=0;i<4;i++)    dx02[i]=f1[i](g_ij,gij_,x0,x1,x2,t+0.5*dt)*dt;
for(int i=0;i<4;i++)   dx12[i]=f1[i](g_ij,gij_,x0,x1,x2,t+0.5*dt)*dt;
for(int i=0;i<4;i++)    dx22[i]=f2[i](g_ij,gij_,x0,x1,x2,t+0.5*dt)*dt;
  
for(int i=0;i<4;i++)   x0[i] += 0.5*dx02[i]; 
for(int i=0;i<4;i++)   x1[i] += 0.5*dx12[i]; 
for(int i=0;i<4;i++)   x2[i] += 0.5*dx22[i]; 
  
for(int i=0;i<4;i++)    dx03[i]=f1[i](g_ij,gij_,x0,x1,x2,t+0.5*dt)*dt;
for(int i=0;i<4;i++)    dx13[i]=f1[i](g_ij,gij_,x0,x1,x2,t+0.5*dt)*dt;
for(int i=0;i<4;i++)    dx23[i]=f2[i](g_ij,gij_,x0,x1,x2,t+0.5*dt)*dt;
  
for(int i=0;i<4;i++)  x0[i] += dx03[i]; 
for(int i=0;i<4;i++)    x1[i] += dx13[i]; 
for(int i=0;i<4;i++)    x2[i] += dx23[i]; 
  
for(int i=0;i<4;i++)    dx04[i]=f1[i](g_ij,gij_,x0,x1,x2,t+dt)*dt;
for(int i=0;i<4;i++)    dx14[i]=f1[i](g_ij,gij_,x0,x1,x2,t+dt)*dt;
for(int i=0;i<4;i++)   dx24[i]=f2[i](g_ij,gij_,x0,x1,x2,t+dt)*dt;
  
for(int i=0;i<4;i++)    x0[i]+=(dx01[i]+2*(dx02[i]+dx03[i])+dx04[i])/6;
for(int i=0;i<4;i++)  x1[i]+=(dx11[i]+2*(dx12[i]+dx13[i])+dx14[i])/6;
for(int i=0;i<4;i++)  x2[i]+=(dx21[i]+2*(dx22[i]+dx23[i])+dx24[i])/6;
  
}
