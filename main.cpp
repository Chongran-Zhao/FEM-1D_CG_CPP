#include <iostream>
#include "Tools.h"
#include <math.h>
#include <Eigen/Eigen>
#define pi acos(-1)
using namespace std;
using namespace Eigen;
int main()
{
  // physical domain
  const double omega_l = 0.0; 
  const double omega_r = 2.0*pi;

  const double g0 = Func_exact(omega_l); //Dirichlet BC
  const double gL = Func_exact_x(omega_r); //Natural BC

  // number of elements
  const int nElem = 10;

  // quadrature rule
  const int nqp = 10;
  double * qp = new double[nqp];
  double * wq = new double[nqp];
  Gauss(nqp, -1.0, 1.0, qp, wq);

  // polynomial (FEM basis function) degree
  int pp = 1;

  const int nLocBas = pp + 1; // number of local basis, local means element

  const int nFunc = pp * nElem + 1;

  int IEN[nLocBas][nElem];
  for (int ee = 0; ee < nElem; ee++)
  {
    for (int aa = 0; aa < nLocBas; aa++)
    {
      IEN[aa][ee] = ee * pp + aa; 
    } 
  }

  const double hh = (omega_r - omega_l) / double(pp*nElem);

  double x_coor[pp*nElem+1];
  for (int ii = 0; ii < pp*nElem+1; ii++)
  {
    x_coor[ii] = omega_l + double(ii) * hh;
  }

  int ID[nFunc];
  for (int ii = 0; ii < nFunc; ii++)
  {
    ID[ii] = ii;
  }
  ID[0] = -1; // assign the ID for the Dirichlet node to be -1

  SparseMatrix <double,RowMajor> K(nFunc,nFunc);
  VectorXd F(nFunc);
  VectorXd u(nFunc);

  for (int ee = 0; ee < nElem; ee++)
  {
    double k_ele[nLocBas][nLocBas];
    double f_ele[nLocBas];
    double x_ele[nLocBas];

    for (int ii = 0; ii < nLocBas; ii++)
    {
      for (int jj = 0; jj < nLocBas; jj++)
      {
        k_ele[ii][jj] = 0.0;
      }
      f_ele[ii] = 0.0;
      x_ele[ii] = 0.0;
    }
    for (int aa = 0; aa < nLocBas; aa++)
    {
      x_ele[aa] = x_coor[IEN[aa][ee]];
    }

    for (int qua = 0; qua < nqp; qua++)
    {
      double dx_dxi = 0.0;
      double x_qua = 0.0;

      for (int aa = 0; aa < nLocBas; aa++)
      {
        x_qua  = x_qua  + x_ele[aa] * PolyBasis(pp, aa+1, 0, qp[qua]);
        dx_dxi = dx_dxi + x_ele[aa] * PolyBasis(pp, aa+1, 1, qp[qua]);
      }

      double dxi_dx = 1.0 / dx_dxi;

      // element assembly for the K and F
      for (int aa = 0; aa < nLocBas; aa++)
      {
        double Na    = PolyBasis(pp, aa+1, 0, qp[qua]);
        double Na_xi = PolyBasis(pp, aa+1, 1, qp[qua]);

        f_ele[aa] = f_ele[aa] + wq[qua] * Func_source(x_qua) * Na * dx_dxi;

        for (int bb = 0; bb < nLocBas; bb++)
        {
          double Nb_xi = PolyBasis(pp, bb+1, 1, qp[qua]);
          k_ele[aa][bb] = k_ele[aa][bb] + wq[qua] * Na_xi * Nb_xi * dxi_dx;
        }
      } 
    }

    for (int aa = 0; aa < nLocBas; aa++)
    {
      int AA = ID[IEN[aa][ee]];
      if(AA >= 0) 
      {
        F(AA) += f_ele[aa];
        for (int bb = 0; bb < nLocBas; bb++)
        {
          int BB = IEN[bb][ee];
          K.coeffRef(AA,BB) += k_ele[aa][bb];
        }
      }
      else
      {
        K.coeffRef(IEN[aa][ee],IEN[aa][ee]) = 1.0;
        F(IEN[aa][ee]) = g0;
      }
    }
  }
  F(nFunc-1) += gL;

  SparseLU<SparseMatrix<double> > solver;
  solver.compute(K);
  u = solver.solve(F);
  cout << "The number of elements = " << nElem << endl;
  cout << "The degree of Legrange-shape function = " << pp << endl;
  cout << "Displacement u = " << endl;
  cout << u << endl;
}

