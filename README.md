# Intro
A CPP code developed to solve 1D pure-diffusion problem in FEM using continuous Galerkin method.

# Phsical problem
$$u_{,xx} + f(x) = 0 \Sigma = (0,L)$$
$$u(0) = u_{exact}(0)$$
$$u_{,x}(L) = u_{exact,x}(L)$$

# Exact solution
$$u_{exact} = sin(x)$$

# Code details
`void Gauss(int N, double a, double * qp, double * wq)` assigns values to `qp` and `wq`.

`double Polybasis(const int &degree, int i, int der, double &x)` returns specific value of Lagrange-shape function.

`double Func_source(const double &x)` returns value of $f(x)$.

`double Func_exact(const double &x)` returns value of exact solution but not used.

`double Func_exact_x(const double &x)` returns value of 1st derivative of exact solution but not used.

# Run the code
1. Run `mkdir build` under `FEM-1D_CG_CPP/` to create a new directory.
2. Run `cd build` to step in the new directory.
3. Run `cmake ..` to create `makefile` using `CMakeLists.txt` from upper directory.
4. Run `make` to create an executable file named `1D_CG_CPP`.
5. Run `./1D_CG_CPP`.
