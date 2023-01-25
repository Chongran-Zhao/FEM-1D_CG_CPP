#include <iostream>
#include "Tools.h"
#include <math.h>
#define pi acos(-1)
#define eps 2.2204e-16
using namespace std;

// u_xx = f
double Func_source(const double & x)
{
  return sin(x);
}

// exacr solution
double Func_exact(const double & x)
{
  return sin(x);
}

// 1st order of derivative of exact solution
double Func_exact_x(const double & x)
{
  return cos(x);
}

void Gauss(int N, double a, double b, double * qp, double * wq)
{
  N -= 1;
  const int N1 = N + 1;
  const int N2 = N + 2;
  double xu[N1];
  for (int ii = 0; ii < N1; ii++)
  {
    xu[ii] = -1.0 + double(ii) * 2.0/(double(N1)-1.0);
  }
  double y[N1];
  for (int ii = 0; ii < N1; ii++)
  {
    y[ii] = cos((2.0*double(ii)+1)*pi/(2.0*double(N)+2.0))+(0.27/double(N1))*sin(pi*xu[ii]*double(N)/double(N2));
  } 
  double L[N1][N2];
  double * Lp = new double[N1];
  double y0[N1];
  double error = 1.0;

  while (error > eps)
  {
    for (int ii = 0; ii < N1; ii++)
    {
      L[ii][0] = 1.0;
      L[ii][1] = y[ii]; 
    }
    for (int ii = 1; ii < N1; ii++)
    {
      for (int jj = 0; jj < N1; jj++)
      {
        L[jj][ii+1] = ( (2.0*double(ii)+1.0) * y[jj] * L[jj][ii] - double(ii)*L[jj][ii-1] ) / double(ii+1);
      }
    }
    for (int ii = 0; ii < N1; ii++)
    {
      Lp[ii] = double(N2) * (L[ii][N1-1] - y[ii]*L[ii][N2-1] ) / (1.0 - y[ii]*y[ii]);
    }
    for (int ii = 0; ii < N1; ii++)
    {
      y0[ii] = y[ii];
      y[ii]  =y0[ii] - L[ii][N2-1] / Lp[ii];
    }

    double error0 = 0.0;
    for (int ii = 0; ii < N1; ii++)
    {
      error = (error0 > abs( y[ii]-y0[ii] )) ? error : abs(y[ii]-y0[ii]);
    }
    error0 = error;
  }
  for (int ii = 0; ii < N1; ii++)
  {
    qp[ii] = (a*(1.0-y[ii])+b*(1.0+y[ii])) / 2.0;
    wq[ii] = (b-a) / ((1-y[ii]*y[ii]) * Lp[ii]*Lp[ii]) * (double(N2)/double(N1)) * (double(N2)/double(N1));

  }
  return;
}


double PolyBasis(const int &degree, int i, int der, double &x)
{
  double poly = 0.0;
  switch (degree)
  {
    // linear basis function
    case 1:
      if (i == 1)
      {
        if (der == 0)
          poly = 0.5 * (1.0-x);
        else if (der == 1)
          poly = -0.5;
        else if (der == 2)
          poly = 0.0;
      }

      if (i == 2)
      {
        if (der == 0)
          poly = 0.5 * (1.0+x);
        else if (der == 1)
          poly = 0.5;
        else if (der == 2)
          poly = 0.0;
      }
      break;

      // quadratic basis function
    case 2:
      if (i == 1)
      {
        if (der == 0)
          poly = 0.5 * x * (x-1.0);
        else if (der == 1)
          poly = x - 0.5;
        else if (der == 2)
          poly = 1.0;
      }

      if (i == 2)
      {
        if (der == 0)
          poly = 1.0 - x*x;
        else if (der == 1)
          poly = -2.0 * x;
        else if (der == 2)
          poly = -2.0;
      }

      if (i == 3)
      {
        if (der == 0)
          poly = 0.5 * x * (x+1.0);
        else if (der == 1)
          poly = x + 0.5;
        else if (der == 2)
          poly = 1.0;
      }

      // Cubic basis function
    case 3:
      if (i == 1)
      {
        if (der == 0)
          poly = -9.0 *(x-(1.0/3.0)) * (x+(1.0/3.0)) * (x-1.0)/16.0;
        else if (der == 1)
          poly = -9.0 * (2.0*x*(x-1.0) + x*x -(1.0/9.0)) / 16.0;
        else if (der == 2)
          poly = -27.0/8.0*x+9.0/8.0;
      }

      if (i == 2)
      {
        if (der == 0)
          poly = 27.0 *(x*x - 1.0)*(x-(1.0/3.0))/16.0;
        else if (der == 1)
          poly = 27.0 * (2.0*x*(x-(1.0/3.0))+x*x-1.0)/16.0;
        else if (der == 2)
          poly = 81.0/8.0*x-9.0/8.0;
      }

      if (i == 3)
      {
        if (der == 0)
          poly = -27.0 * (x*x-1.0)*(x+(1.0/3.0))/16.0;
        else if (der == 1)
          poly = -27.0 * (2.0*x*(x+(1.0/3.0))+x*x-1.0)/16.0;
        else if (der == 2)
          poly = -81.0/8.0*x-9.0/8.0;
      }

      if (i == 4)
      {
        if (der == 0)
          poly = 9.0*(x+1.0)*(x*x-(1.0/9.0))/16.0;
        else if (der == 1)
          poly = 9.0*(x*x-(1.0/9.0)+2.0*x*(x+1.0))/16.0;
        else if (der == 2)
          poly = 27.0/8.0*x+9.0/8.0;
      }

      // quartic basis function
    case 4:
      if (i == 1)
      {
        if (der == 0)
          poly = 2.0*x*(x*x-(1.0/4.0))*(x-1.0)/3.0;
        else if (der == 1)
          poly = 2.0*((x*x-(1.0/4.0))*(x-1.0)+2.0*x*x*(x-1.0)+x*(x*x-(1.0/4.0)))/3.0;
        else if (der == 2)
          poly = 4.0*x*(x-1)+4.0*x*x-1.0/3.0;
      }

      if (i == 2)
      {
        if (der == 0)
          poly = -8.0*x*(x*x-1.0)*(x-0.5)/3.0;
        else if (der == 1)
          poly = -8.0*((x*x-1.0)*(x-0.5)+x*x*(2.0*x-1.0)+x*(x*x-1.0))/3.0;
        else if (der == 2)
          poly = -16.0/3.0*x*(x-0.5)-16.0*x*x+16.0/3.0-16.0/3.0*x*(2.0*x-1.0);
      }

      if (i == 3)
      {
        if (der == 0)
          poly = 4.0*(x*x-1.0)*(x*x-0.25);
        else if (der == 1)
          poly = 4.0*(2.0*x*(x*x-0.25)+2.0*x*(x*x-1.0));
        else if (der == 2)
          poly = 48.0*x*x-10.0;
      }

      if (i == 4)
      {
        if (der == 0)
          poly = -8.0*x*(x*x-1.0)*(x+0.5)/3.0;
        else if (der == 1)
          poly = -8.0*((x*x-1.0)*(x+0.5)+x*x*(2.0*x+1.0)+x*(x*x-1.0))/3.0;
        else if (der == 2)
          poly = -16.0/3.0*x*(x+0.5)-16.0*x*x+16.0/3.0-16.0/3.0*x*(2.0*x+1.0);
      }

      if (i == 5)
      {
        if (der == 0)
          poly = 2.0*x*(x*x-0.25)*(x+1.0)/3.0;
        else if (der == 1)
          poly = 2.0*((x*x-0.25)*(x+1.0)+2.0*x*x*(x+1.0)+x*(x*x-0.25))/3.0;
        else if (der == 2)
          poly = 4.0*x*(x+1.0)+4.0*x*x-1.0/3.0;
      }

      // quintic basis function
    case 5:
      if (i == 1)
      {
        if (der == 0)
          poly =-625.0*(x*x-(9.0/25.0))*(x*x-(1.0/25.0))*(x-1.0)/768.0;
        else if (der == 1)
          poly = -3125.0/768.0*x*x*x*x+625.0/192.0*x*x*x+125.0/128.0*x*x-125.0/192.0*x-3.0/256.0;
        else if (der == 2)
          poly = -3125.0/192.0*x*x*x+625.0/64.0*x*x+125.0/64.0*x-125.0/192.0;
      }

      if (i == 2)
      {
        if (der == 0)
          poly = 3125.0/768.0*(x+1.0)*(x+0.2)*(x-0.2)*(x-0.6)*(x-1.0);
        else if (der == 1)
          poly = 15625.0/768.0*x*x*x*x-625.0/64.0*x*x*x-1625.0/128.0*x*x+325.0/64.0*x+125.0/768.0;
        else if (der == 2)
          poly = 15625.0/192.0*x*x*x-1875.0/64.0*x*x-1625.0/64.0*x+325.0/64.0;
      }

      if (i == 3)
      {
        if (der == 0)
          poly = -3125.0/384.0*(x+1.0)*(x+0.6)*(x-0.2)*(x-0.6)*(x-1.0);
        else if (der == 1)
          poly = -15625.0/384.0*x*x*x*x+625.0/96.0*x*x*x+2125.0/64.0*x*x-425.0/96.0*x-375.0/128.0;
        else if (der == 2)
          poly = -15625.0/96.0*x*x*x*x+625.0/32.0*x*x+2125.0/32.0*x-425.0/96.0;
      }

      if (i == 4)
      {
        if (der == 0)
          poly = 3125.0/384.0*(x+1.0)*(x+0.6)*(x+0.2)*(x-0.6)*(x-1.0);
        else if (der == 1)
          poly = 15625.0/384.0*x*x*x*x+625.0/96.0*x*x*x-2125.0/64.0*x*x-425.0/96.0*x+375.0/128.0;
        else if (der == 2)
          poly = 15625.0/96.0*x*x*x+625.0/32.0*x*x-2125.0/32.0*x-425.0/96.0;
      }

      if (i == 5)
      {
        if (der == 0)
          poly = -3125.0/768.0*(x+1.0)*(x+0.6)*(x+0.2)*(x-0.2)*(x-1.0);
        else if (der == 1)
          poly = -15625.0/768.0*x*x*x*x-625.0/64.0*x*x*x+1625.0/128.0*x*x+325.0/64.0*x-125.0/768.0;
        else if (der == 2)
          poly = -15625.0/192.0*x*x*x-1875.0/64.0*x*x+1625.0/64.0*x+325.0/64.0;
      }

      if (i == 6)
      {
        if (der == 0)
          poly = 625.0/768.0*(x+1.0)*(x+0.6)*(x+0.2)*(x-0.2)*(x-0.6);
        else if (der == 1)
          poly = 3125.0/768.0*x*x*x*x-125.0/128.0*x*x+3.0/256.0+625.0/192.0*x*x*x-125.0/192.0*x;
        else if (der == 2)
          poly = 3125.0/192.0*x*x*x-125.0/64.0*x+625.0/64.0*x*x-125.0/192.0;
      }

    case 6:
      if (i == 1)
      {
        if (der == 0)
          poly = 81.0/80.0*(x+2.0/3.0)*(x+1.0/3.0)*x*(x-1.0/3.0)*(x-2.0/3.0)*(x-1.0);
        else if (der == 1)
          poly = 1.0/80.0*(6.0*x-1.0)*(81.0*x*x*x*x-54.0*x*x*x-39.0*x*x+16.0*x+4.0);
        else if (der == 2)
          poly = 243.0/40.0*x*x*x*x-81.0/20.0*x*x*x-117.0/40.0*x*x+6.0/5.0*x+3.0/10.0+(3.0/40.0*x-1.0/80.0)*(324.0*x*x*x-162.0*x*x-78.0*x+16.0);
      } 

      if (i == 2)
      {
        if (der == 0)
          poly = -243.0/40.0*(x+1.0)*(x+1.0/3.0)*x*(x-1.0/3.0)*(x-2.0/3.0)*(x-1.0);
        else if (der == 1)
          poly = -27.0/20.0*x-27.0/2.0*x*x+81.0/4.0*x*x*x*x-729.0/20.0*x*x*x*x*x+27.0*x*x*x+9.0/20.0;
        else if (der == 2)
          poly = -27.0/20.0-27.0*x+81.0*x*x*x-729.0/4.0*x*x*x*x+81.0*x*x;
      }

      if (i == 3)
      {
        if (der == 0)
          poly = 243.0/16.0*(x+1.0)*(x+2.0/3.0)*x*(x-1.0/3.0)*(x-2.0/3.0)*(x-1.0);
        else if (der == 1)
          poly = 27.0/2.0*x+351.0/16.0*x*x-405.0/16.0*x*x*x*x-2.25+729.0/8.0*x*x*x*x*x-351.0/4.0*x*x*x;
        else if (der == 2)
          poly = 27.0/2.0+351.0/8.0*x-405.0/4.0*x*x*x+3645.0/8.0*x*x*x*x-1053.0/4.0*x*x;
      }

      if (i == 4)
      {
        if (der == 0)
          poly = -81.0/4.0*(x+1.0)*(x+1.0/3.0)*(x+2.0/3.0)*(x-1.0/3.0)*(x-2.0/3.0)*(x-1.0);
        else if (der == 1)
          poly = -49.0/2.0*x-243.0/2.0*x*x*x*x*x+126.0*x*x*x;
        else if (der == 2)
          poly = -49.0/2.0-1215.0/2.0*x*x*x*x+378.0*x*x;
      }

      if (i == 5)
      {
        if (der == 0)
          poly = 243.0/16.0*(x+1.0)*(x+1.0/3.0)*(x+2.0/3.0)*x*(x-2.0/3.0)*(x-1.0);
        else if (der == 1)
          poly = 27.0/2.0*x-351.0/16.0*x*x+405.0/16.0*x*x*x*x+729.0/8.0*x*x*x*x*x+2.25-351.0/4.0*x*x*x;
        else if (der == 2)
          poly = 27.0/2.0-351.0/8.0*x+405.0/4.0*x*x*x+3645.0/8.0*x*x*x*x-1053.0/4.0*x*x;
      }

      if (i == 6)
      {
        if (der == 0)
          poly = -243.0/40.0*(x+1.0)*(x+1.0/3.0)*(x+2.0/3.0)*x*(x-1.0/3.0)*(x-1.0);
        else if (der == 1)
          poly = -27.0/20.0*x+27.0/2.0*x*x-81.0/4.0*x*x*x*x-729.0/20.0*x*x*x*x*x+27.0*x*x*x-9.0/20.0;
        else if (der == 2)
          poly = -27.0/20.0+27.0*x-81.0*x*x*x-729.0/4.0*x*x*x*x+81.0*x*x;
      }

      if (i == 7)
      {
        if (der == 0)
          poly = 81.0/80.0*(x+1.0)*(x+1.0/3.0)*(x+2.0/3.0)*x*(x-1.0/3.0)*(x-2.0/3.0);
        else if (der == 1)
          poly = 1.0/80.0*(6.0*x+1.0)*(81.0*x*x*x*x+54.0*x*x*x-39.0*x*x-16.0*x+4.0);
        else if (der == 2)
          poly = 243.0/40.0*x*x*x*x+81.0/20.0*x*x*x-117.0/40.0*x*x-6.0/5.0*x+3.0/10.0+(3.0/40.0*x+1.0/80.0)*(324.0*x*x*x+162.0*x*x-78.0*x-16.0);
      }    
  }
  return poly;
}


