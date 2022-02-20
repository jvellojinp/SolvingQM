#include <iostream>
#include <vector>
#include <algorithm>
#include <cstdlib>
#include <math.h>


double rk4(double ta, double tb, double h, std::vector<double> & y, double W);
double f(double t, const std::vector<double> & y, int id, double W);

int main (void)
{
  const double N  = 2;
  const double TA = 0.001;
  const double TB = 0.999;
  const double H  = 0.01;
  const double HW = 0.001; 
  double W = 0.0; // omega, in rad/s
  std::vector<double> y = {0, 1.0}; // {x0, v0}

  for (int i = 0; i < 20000; ++i)
    {
      std::cout << W << '\t' << rk4(TA,TB,H,y,W) << '\n';

      W = W + HW;
      y[0]=0;
      y[1]=1.0;
    }

  return 0;
}

double rk4(double ta, double tb, double h, std::vector<double> & y,double W)
{
  std::vector<double> k1, k2, k3, k4, aux;
  k1.resize(y.size());
  k2.resize(y.size());
  k3.resize(y.size());
  k4.resize(y.size());
  aux.resize(y.size());

  const int N = (tb-ta)/h;
  for (int nt = 0; nt < N; ++nt) {
    double t = ta + h*nt;
    // k1
    for(int ii = 0; ii < y.size(); ++ii) {
      k1[ii] = h*f(t, y, ii, W);
    }
    // k2 aux
    for(int ii = 0; ii < y.size(); ++ii) {
      aux[ii] = y[ii] + k1[ii]/2;
    }
    //k2
    for(int ii = 0; ii < y.size(); ++ii) {
      k2[ii] = h*f(t + h/2, aux, ii, W);
    }
    // k3 aux
    for(int ii = 0; ii < y.size(); ++ii) {
      aux[ii] = y[ii] + k2[ii]/2;
    }
    //k3
    for(int ii = 0; ii < y.size(); ++ii) {
      k3[ii] = h*f(t + h/2, aux, ii, W);
    }
    // k4 aux
    for(int ii = 0; ii < y.size(); ++ii) {
      aux[ii] = y[ii] + k3[ii];
    }
    //k4
    for(int ii = 0; ii < y.size(); ++ii) {
      k4[ii] = h*f(t + h, aux, ii, W);
    }
    // write new y
    for(int ii = 0; ii < y.size(); ++ii) {
      y[ii] = y[ii] + (k1[ii] + 2*k2[ii] + 2*k3[ii] + k4[ii])/6.0;
    }
    //std::cout << t << "\t" << y[0] << "\t" << y[1] << std::endl;
  }
  return y[0];
}


double f(double t, const std::vector<double> & y, int id, double W)
{
  if (0 == id) {
    return y[1];
  }
  else if (1 == id) {
    return (100*std::exp(-100*(1.0-2.0*t)*(1.0-2.0*t))-W*W)*y[0];// - 0.23*y[1];
  }
  else {
    std::cerr << "ERROR!!!!!!!! Id does not exists -> " << id << std::endl;
    exit(1);
  }
}
