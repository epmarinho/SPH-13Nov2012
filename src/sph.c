/**@<sph.c>::*/
#include <string.h>
#include <math.h>

#include <constants.h>
#include <k-nn.h>
#include <treesph.h>
#include <kernel.h>
#include <sph.h>
#include <macros.h>

#define Rgas    8.314e+07 /* erg mole⁻¹ Kelvin⁻¹ */
#define adbtc   .4 /* valid for diatomic gases only */
#define Av      6.0221415e+23 /* particles mole⁻¹ */
#define H       1 / 6.0221415e+23/* mass in grams of 1 mole of H atoms */
#define H2      2 / 6.0221415e+23/* mass in grams of 1 mole of molecular H */
#define O2      32 / 6.0221415e+23/* mass in grams of 1 mole of molecular O */
#define rhoH2   2 / 22.414e+03 /* g cm⁻³ for H2 */
#define rhoO2   32 / 22.414e+03 /* g cm⁻³ for O2 */
#define u0      UMIN

double
e_crit          = .5,
rho_crit        = 10,
rho_1           = .01,
rho_2           = .1,
rho_3           = 1,
Gamma           = 2,
Lambda_1        = 1,
Lambda_2        = 10;

extern double Etot, Kinect, Ugrav, Utherm;
extern double Phi[NMAX];

double          rho[NMAX];
void
sph_densities(int N)
{
  int             i;
  init_kernels(N);
  for (i = 0; i < N; i++) {
    int             l;
    for (rho[i] = l = 0; l < KNN[i].population; l++) {
      int             j = KNN[i].list[l] - 1;
      rho[i] += mass[j] * W[i][l];
    }
  }
}

double
u[NMAX],
udot[NMAX],
Paccel[NMAX][DIM];

/*
 * artificial viscosity parameters (Monaghan-Gingold 1983; Marinho-Lépine
 * 2000)
 */
double
alpha   = 8,
beta    = 20,
eta     = .0625;

/*
 * Most SPH quantities are estimated here
 */
void
sph_quantities(int n)
{
  int             i;

  sph_densities(n);

  #ifndef _INTEGRATE_ON_ENERGIES_
  for (i = 0; i < n; i++){
    u[i] = max(min(u0 * pow(rho[i] / rhoH2, adbtc), UMAX), UMIN);
  }
  #endif

  for (i = 0; i < n; i++) {

    bzero(Paccel[i], sizeof(double) * DIM);

    udot[i] = 0;

    int             l;
    for (l = 1; l < KNN[i].population; l++) {

      int             j = KNN[i].list[l] - 1;
      int k;

      /*
       * symmetrizing SPH-pressure tensor
       */
      double          Pi_ij;

      /*
       * symmetrizing SPH-pressure model
       */
      Pi_ij = adbtc * (u[j] / rho[j] + u[i] / rho[i]);

      /*
       * artificial viscosity model
       */
      double
      mu_ij,
      h_ij = .5 * (h[i] + h[j]),
      r_ij = sqrt(query2point_distance(x[i], x[j]));
      for (mu_ij = k = 0; k < DIM; k++) {
        mu_ij += (v[i][k] - v[j][k]) * (x[i][k] - x[j][k]);
      }
      if (mu_ij < 0) {
        mu_ij =
        -mu_ij / (r_ij * (r_ij / h_ij) + eta * h_ij);
        Pi_ij +=
        (alpha * sqrt(u[i] + u[j]) +
        beta * mu_ij) * mu_ij / (rho[i] + rho[j]);
      }
      //end artificial viscosity model

      /*
       * energy evolution equation
       */
      double          B_jij = mass[j] * Pi_ij;
      for (k = 0; k < DIM; k++) {
        double          C_ijk = gradW[i][l][k] * B_jij;

        /*
         * get ride of the thermal energy computation to get also
         * pressure acceleration:
         */
        Paccel[i][k] += C_ijk;

        /*
         * adiabatic heating:
         */
        udot[i] += .5 * (v[i][k] - v[j][k]) * C_ijk;

      }//endfor k
    }//endfor l

#ifdef __ADOPT_CRITICAL_COUNTERMEASURES__
    /*
     * artificial heating:
     */
    if (e_crit < u[i] && u[i] < UMAX && rho[i] > rho_crit ) {
      udot[i] += Gamma;// * rho[i] / Mtot;
    }

    /*
     * artificial cooling:
     */
    if (u[i] > UMIN) {
      if (rho[i] < rho_1)
        udot[i] -= Lambda_1 * u[i] * Mtot / rho[i];
      else if (rho_2 < rho[i] && rho[i] < rho_3)
        udot[i] -= Lambda_2 * u[i];
    }
#endif

  Utherm += u[i]*mass[i];

  }//endfor i
}







/** The following is used for a naïve way of computing SPH quantities, as it was
 *    first introduced in both original papers of Gingold and Monaghan (1977)
 *    and of Lucy (1977) **/


#ifdef _FOR_TESTS_ONLY_
double          Bias[NMAX][DIM];
void
sph_bias_field(int n)
{
  int             i;
  for (i = 0; i < n; i++) {
    int             l,
    k;
    bzero(Bias[i], sizeof(double) * DIM);
    for (l = 1; l < KNN[i].population; l++) {
      int             j = KNN[i].list[l] - 1;
      for (k = 0; k < DIM; k++) {
        Bias[i][k] +=
        gradW[i][l][k] * mass[j] * (1 / pow(rho[j], 2) +
        1 / pow(rho[i], 2));
      }
    }
    for (k = 0; k < DIM; k++) {
      Bias[i][k] *= (rho[i] / norm[i]);
    }
  }
}
double          rho_gradient[NMAX][DIM];
void
sph_density_gradient(int n)
{
  int             i;
  for (i = 0; i < n; i++) {
    int             l,
    k;
    bzero(rho_gradient[i], sizeof(double) * DIM);
    for (l = 1; l < KNN[i].population; l++) {
      int             j = KNN[i].list[l] - 1;
      for (k = 0; k < DIM; k++) {
        rho_gradient[i][k] +=
        gradW[i][l][k] * mass[j] * (1 / rho[j] + 1 / rho[i]);
      }
    }
    for (k = 0; k < DIM; k++) {
      rho_gradient[i][k] /= norm[i];
      rho_gradient[i][k] -= Bias[i][k];
    }
  }
}

double
smoothing(int i, double A[])
{
  double          Asmooth = 0;
  int             l;
  for (l = 0; l < KNN[i].population; l++) {
    int             j = KNN[i].list[l] - 1;
    Asmooth += W[i][l] * A[j] * mass[j] / rho[j];
  }
  return Asmooth;
}

double
smooth_derivative(int i, int k, double A[])
{
  double          gradA = 0;
  int             l;
  for (l = 0; l < KNN[i].population; l++) {
    int             j = KNN[i].list[l] - 1;
    gradA += gradW[i][l][k] * A[j] * mass[j] / rho[j];
  }
  return gradA;
}

double          norm[NMAX];
#endif
