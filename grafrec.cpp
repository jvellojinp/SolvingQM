#include <iostream>
#include <vector>
#include <algorithm>
#include <cstdlib>
#include <math.h>
#include <fstream>

#define EPSILON 0.0001


using namespace std;


double f1(double t, double x, double v)
{
  return(v);
}

double f2(double t, double x, double v, double k)
{
  return( (100*pow(2.7,-(1-2*x)*(1-2*x))-k*k)*x );
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

// Bisection function
double bisection(double a, double b)
{
  double t0 = 0, x0 = 0, v0 = 1, h = 0.001;
  if (r(t0, x0, v0, 1, h, a)*r(t0, x0, v0, 1, h, b) >= 0)
  {
      cout << "You have not assumed right a and b\n";
  }

  double c;
  while ((b-a) >= EPSILON)
  {
      // Find middle point
      c = (a+b)/2;

      // Check if middle point is root
      if (r(t0, x0, v0, 1, h, c) == 0.0)
          break;

      // Decide the side to repeat the steps
      else if (r(t0, x0, v0, 1, h, c)*r(t0, x0, v0, 1, h, a) < 0)
          b = c;
      else
          a = c;
  }
  return c;
}

int main()
{
  int root_counter = 0, root_number = 5;
  double k = 0, k_next = 0, t0 = 0, x0 = 0, v0 = 1, h = 0.001, step = 0.001;

  do
  {
      //Inittialy try with step = .01 and root_counter = 1,2
      
      k_next = k + step;
      if(r(t0, x0, v0, 1, h, k)*r(t0, x0, v0, 1, h, k_next) <= 0)
      {
          root_counter = root_counter + 1;
          cout << "Root interval #" << root_counter << ": " << "a = " << k << ", b = " << k_next << "-> Root = " << bisection(k,k_next) <<  endl;
      }
      k = k_next;
  }
  while(k <= 20  && root_counter != root_number);

  return 0;
}
