#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <strings.h>
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
double          Damp = 1;
void
octree_gravity_rotation(int n)
{
    int             i,
                    k;

    for (i = 0; i < n; i++) {
        int             k;
        double          r,
                        vc;
    /** perform the topdown search for the nearest neighbors seen inside
     * the solid angle Omega */
        Octree_descent(i);

        for (bzero(g[i], sizeof(double) * DIM), k = 0; k < accepted; k++) {
            int             nu = welldistant[k],
                j;
            double          r =
                distance(Octree[nu].vertex, x[i]) + epsilon;

            for (j = 0; j < DIM; j++) {
                double          delta = Octree[nu].vertex[j] - x[i][j];
                g[i][j] += Octree[nu].mass * delta / (r * sqrt(r));
            }

        }

        r = sqrt(x[i][0] * x[i][0] + x[i][1] * x[i][1]);
        vc = sqrt(g[i][0] * g[i][0] + g[i][1] * g[i][1]);
        vc = Damp * sqrt(vc / r);

        v[i][0] += vc * (-x[i][1]);
        v[i][1] += vc * (x[i][0]);

    }
}

main(int argc, char const *argv[])
{
    int             n = NMAX,
        i;
    FILE           *input;

    if (argc > 1) {
        n = atoi(argv[1]);
        input = fopen(argv[1], "r");
    } else {
        input = stdin;
    }

    if (argc > 2) {
        Damp = atof(argv[2]);
    }

    if (argc > 3) {
        Omega = atof(argv[3]);
    }

    if (argc > 4) {
        epsilon = atof(argv[4]);
    } else {
        epsilon = .01;
    }

    Octree_build(n = getdata(n, input));

    octree_gravity_rotation(n);

    show_result(n, stdout);

    exit(0);
}
