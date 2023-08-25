/**@<sph.h>::*/
#pragma once

#define UMAX    10000
#define UMIN    .001

extern double u[NMAX], udot[NMAX], P[NMAX], Paccel[NMAX][DIM];
extern void sph_quantities(int n);
extern double rho[NMAX];
extern void sph_densities(int N);
extern double smoothing(int i, double A[]);
extern double smooth_derivative(int i, int k, double A[]);
extern void sph_bias_field(int n);
extern void sph_density_gradient(int n);
extern double Bias[NMAX][DIM];
extern double rho_gradient[NMAX][DIM];
