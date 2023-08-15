/**@<leapfrog.c>::**/
#include <stdlib.h>
#include <stdio.h>
#include <strings.h>
#include <math.h>
#include <ctype.h>

#include <macros.h>
#include <constants.h>
#include <octree.h>
#include <n-body.h>
#include <io-stuff.h>
#include <k-nn.h>
#include <kernel.h>
#include <sph.h>
#include <leapfrog.h>
double
scalarprod(double x[], double y[])
{
  int             j;
  double          s = 0;
  for (j = 0; j < DIM; j++)
    s += x[j] * y[j];
  return s;
}

double          dtold[NMAX];
leapfrog(int n, double DT)
{
  int             i,
  j;
  double          dtmin = (double) 1 / 0x10000;

  for (i = 0; i < n; i++) {
    double          dt = DT,
    Dt,
    tau,
    dtau,
    vscalar,
    paccell;
    for (Dt = 0; Dt < DT; Dt += dt) {

      /*
       * guess a typical length scale for current particle
       */
      double          ds =
      min(h[i] / pow(KNN[i].population, .3333333333333333),
          sqrt(Octree[crossreference[i]].size));

      /*
       * present length scale should not be crossed by the particle
       * along a reasonable time scale
       */
      while (vscalar = sqrt(scalarprod(v[i], v[i])) * dt > ds
        && dt > dtmin)
        dt /= 2;

      /*
       * particle accelleration cannot gain velocity enough to cross
       * the velocity scale
       */
      while (vscalar
        && (paccell =
        sqrt(scalarprod(Paccel[i], Paccel[i])) * dt) >
        vscalar && dt > dtmin)
        dt /= 2;

      /*
       * similar to thermal energy rate
       */
      if (udot[i] != 0) {
        while (fabs(u[i] / udot[i]) < dt && dt > dtmin)
          dt /= 2;
      }

      /*
       * leapfrog here
       */
      tau = (dt + dtold[i]) * .5;
      dtau = (dt - dtold[i]) * .5;
      for (j = 0; j < DIM; j++) {
        double accell = g[i][j] - Paccel[i][j];
        x[i][j] += (v[i][j] + .5 * accell * dtau) * tau;
        v[i][j] += accell * tau;
      }

      #ifdef _INTEGRATE_ON_ENERGIES_
      u[i] += udot[i] * tau;
      /*
       * avoid unreal negative thermal energy and strong explosions
       */
      u[i] = min(max(u[i], UMIN), UMAX);
      #endif
      dtold[i] = dt;
      /*
       * end leapfrog
       */
    }
  }
}
