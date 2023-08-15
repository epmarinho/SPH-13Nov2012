/**@<treesph.c>::*/
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

double          Mtot = 0;
double          mass[NMAX];
double          x[NMAX][DIM];
double          v[NMAX][DIM];
double          g[NMAX][DIM];
double          dt = 1. / 128;  // it depends on the gravity time scale

main(int argc, char const *argv[])
{
  int             n = NMAX,
  i;
  FILE           *input,
  *output;

  if (argc > 1) {
    input = fopen(argv[1], "r");
  } else {
    input = stdin;
  }

  if (argc > 2) {
    n = atoi(argv[2]);
    getdata(n, input);
    Octree_build(n);
  } else {
    Octree_build(n = getdata(n, input));
  }

  if (argc > 3) {
    K = atoi(argv[3]);
  } else {
    K = ceil(sqrt((double) n));
  }
  K = min(K, MAXNEIGHBORINGS);

  if (argc > 4) {
    Omega = atof(argv[4]);
  }

  if (argc > 5) {
    epsilon = atof(argv[5]);
  }

  {
    double          dt2;
    if (argc > 6) {
      dt = atof(argv[6]);
    }
    for (dt2 = 1; dt2 > dt; dt2 *= .5);
    dt = dt2;
  }

  /**
   * first prediction:
   * ##############################################################
   *  @-> gravity and sph are computed from the initial conditions:
   *          <<x(n-1/2) v(n) u(n-1/2)>>
   * ##############################################################
   */

  octree_gravity(n);

  sph_quantities(n);

  leapfrog(n, .5 * dt);

  /**
   * second prediction
   * ######################################################
   *  @-> only sph is computed from the initial conditions:
   *          <<x(n) v(n+1/2) u(n)>>
   * ######################################################
   */

  sph_quantities(n);

  leapfrog(n, .5 * dt);

  /**
   * output
   * ######################################################
   *  @-> final conditions:
   *          <<x(n+1/2) v(n+1) u(n+1/2)>>
   * ######################################################
   */
  show_result(n, output = stdout);

  exit(0);
}
