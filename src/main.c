/**@<main.c>::*/
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
double          Etot, Kinect, Ugrav, Utherm;

extern double Phi[NMAX];


int main(int argc, char const *argv[])
{
  int             n = NMAX,
  i;
  FILE           *input,
  *output, *Energies;

  if (argc > 1) {
    input = fopen(argv[1], "r");
    fprintf(stderr,"input file: %s\n", argv[1]);
  } else {
    input = stdin;
  }

  if (argc > 2) {
    n = atoi(argv[2]);
    fprintf(stderr,"number of particles: %d\n", n);
    fprintf(stderr, "entering getdata\n");
    getdata(n, input);
    fprintf(stderr, "entering Octree_build for n = %d\n", n);
    Octree_build(n);
  } else {
    fprintf(stderr, "entering Octree_build for n = %d\n", n);
    Octree_build(n = getdata(n, input));
  }

  if (argc > 3) {
    K = atoi(argv[3]);
  } else {
    K = ceil(sqrt((double) n));
  }
  K = min(K, MAXNEIGHBORINGS);
  fprintf(stderr,"K = %d\n", K);

  if (argc > 4) {
    Omega = atof(argv[4]);
    fprintf(stderr,"Omega = %lg\n", Omega);
  }

  if (argc > 5) {
    epsilon = atof(argv[5]);
    fprintf(stderr,"epsilon = %lg\n", epsilon);
  }

  {
    double          dt2;
    if (argc > 6) {
      dt = atof(argv[6]);
      fprintf(stderr,"dt = %lg\n", dt);
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

  fprintf(stderr, "entering octree_gravity\n");
  octree_gravity(n);

  fprintf(stderr, "entering sph_quantities\n");
  sph_quantities(n);

  fprintf(stderr, "entering leapfrog\n");
  leapfrog(n, .5 * dt);

  /**
   * second prediction
   * ######################################################
   *  @-> only sph is computed from the initial conditions:
   *          <<x(n) v(n+1/2) u(n)>>
   * ######################################################
   */

  fprintf(stderr, "entering sph_quantities again\n");
  sph_quantities(n);

  fprintf(stderr, "entering leapfrog too\n");
  leapfrog(n, .5 * dt);

  Etot = Kinect + Ugrav + Utherm;

  /**
   * output
   * ######################################################
   *  @-> final conditions:
   *          <<x(n+1/2) v(n+1) u(n+1/2)>>
   * ######################################################
   */
  fprintf(stderr, "entering show_result\n");
  show_result(n, output = stdout);

  Energies = fopen("Enengies_hist.txt", "a+");
  fprintf(Energies,"%f %f %f %f\n",Etot, Kinect, Ugrav, Utherm);
  exit(0);
}
