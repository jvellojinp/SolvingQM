#include <iostream>
#include <vector>
#include <algorithm>
#include <cstdlib>
#include <math.h>
#include <fstream>

using namespace std;

double f1(double t, double x, double v)
{
  return(v);
}

double f2(double t, double x, double v, double k)
{
  return(( 100*x-k*k)*v);
}

double r(double t0, double x0, double v0, double t, double h, double k)
{ 
  int n = (int)((t-t0)/h);

  double k1x, k1v, k2x, k2v, k3x, k3v, k4x, k4v;

  for (int i=1; i<=n; i++)
  {
    k1x = f1(t0, x0, v0);
    k1v = f2(t0, x0, v0, k);
    k2x = f1(t0 + 0.5*h, x0 + 0.5*h*k1x, v0 + 0.5*h*k1v);
    k2v = f2(t0 + 0.5*h, x0 + 0.5*h*k1x, v0 + 0.5*h*k1v, k);
    k3x = f1(t0 + 0.5*h, x0 + 0.5*h*k2x, v0 + 0.5*h*k2v);
    k3v = f2(t0 + 0.5*h, x0 + 0.5*h*k2x, v0 + 0.5*h*k2v, k);
    k4x = f1(t0 + h, x0 + h*k3x, v0 + h*k3v);
    k4v = f2(t0 + h, x0 + h*k3x, v0 + h*k3v, k);

    v0 = v0 + h*(1.0/6.0)*(k1v + 2*k2v + 2*k3v + k4v);
    x0 = x0 + h*(1.0/6.0)*(k1x + 2*k2x + 2*k3x + k4x);
    t0 = t0 + h;
  }

  return x0;
}

int main()
{
  ofstream file;
  file.open("Disparo.txt");
    double t0 = 0, x0 = 0, v0 = 1, h = 0.001, eps = 0.001, k = 0;

    do
    {
      file << k << " " << r(t0, x0, v0, 1, h, k) << endl;
      k = k + 0.001;
    }
    while(k<=20);
  file.close();

  return 0;
}
