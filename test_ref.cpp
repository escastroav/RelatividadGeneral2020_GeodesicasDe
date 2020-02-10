#include <iostream>

using namespace std;

double foo(double (&t)[2])
{
  t[0]+=0.1;
  return 2*t[0];
}

int main()
{
  double t[2] = {0,0};
  
  for(double tt=0;tt<1; ){
    tt = t[0];
    cout << tt << "\t" << foo(t) <<endl;
  }
  
}
