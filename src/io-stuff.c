/**@<io-stuff.c>*/
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <string.h>
#include <math.h>

#include <constants.h>
#include <macros.h>
#include <k-nn.h>
#include <octree.h>
#include <treesph.h>
#include <kernel.h>
#include <sph.h>
#include <leapfrog.h>

skipspaces(FILE * input)
{
    int             c;
    while (isspace(c = getc(input)));
  /** check eof **/
    if (c == -1)
        return -1;
    ungetc(c, input);
    return 0;
}

getdata(int n, FILE * input)
{
    int             j,
                    i;
    for (Mtot = i = 0; 1; i++) {
        int             c;

        skipspaces(input);
        c = fscanf(input, "%le", &mass[i]);
        Mtot += mass[i];

        for (j = 0; j < DIM; j++) {
            skipspaces(input);
            c = fscanf(input, "%le", &x[i][j]);
        }

        for (j = 0; j < DIM; j++) {
            skipspaces(input);
            c = fscanf(input, "%le", &v[i][j]);
        }

        skipspaces(input);
        c = fscanf(input, "%le", &u[i]);

        skipspaces(input);
        c = fscanf(input, "%le", &dtold[i]);

    /** seek eol **/
        while ((c = getc(input)) != '\n');

        if (skipspaces(input))
            return i + 1;
    }
}

double          zero[] = { 0, 0, 0 }, sysradius = 0, radialprojection, r;

void
show_result(int n, FILE * output)
{
    int             i;

    for (i = 0; i < n; i++) {

        int             j;
        double          r;

        // sysradius = max(sysradius, sqrt(query2point_distance(x[i],
        // zero)));

        // 1
        fprintf(output, "%22.14le", mass[i]);

        // 2-4
        for (j = 0; j < DIM; j++)
            fprintf(output, "%22.14le", x[i][j]);

        // 5-7
        for (j = 0; j < DIM; j++)
            fprintf(output, "%22.14le", v[i][j]);

        // 8
        fprintf(output, "%22.14le", u[i]);

        // 9
        fprintf(output, "%22.14le", dtold[i]);

        // 10
        fprintf(output, "%22.14le", udot[i]);

        // 11
        fprintf(output, "%22.14le", rho[i]);

#ifdef _FOR_ERROR_CHECK_ONLY_

        // for (j = 0; j < DIM; j++)
        // fprintf(output, "%22.14le", -Paccel[i][j]);

        // 10-12
        r = sqrt(query2point_distance(x[i], Octree[ROOT].vertex));
        for (radialprojection = j = 0; j < DIM; j++) {
            fprintf(output, "%22.14le", Bias[i][j]);
            radialprojection +=
                Bias[i][j] * (x[i][j] - Octree[ROOT].vertex[j]);
        }

        // 13
        fprintf(output, "%22.14le", radialprojection /= r);

        // 14-16
        for (radialprojection = j = 0; j < DIM; j++) {
            fprintf(output, "%22.14le", rho_gradient[i][j]);
            radialprojection +=
                rho_gradient[i][j] * (x[i][j] - Octree[ROOT].vertex[j]);
        }

        // 17
        fprintf(output, "%22.14le", radialprojection /= r);
#endif
        fprintf(output, "\n");
    }

    // fprintf(stderr, "system radius = %lg\n", sysradius);
}
