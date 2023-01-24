#ifndef _TOOL_H_
#define _TOOL_H_
void Gauss(int N, double a, double b, double * qp, double * wq);

double PolyBasis(const int &degree, int i, int der, double &x);

double Func_source(const double & x);

double Func_exact(const double & x);

double Func_exact_x(const double & x);

#endif
