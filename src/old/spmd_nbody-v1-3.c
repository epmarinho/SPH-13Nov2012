#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <strings.h>
#include <ctype.h>
#include <spmd_octree.h>

double          Mtot = 0;
double          mass[N];
double          x[N][DIM];
double          v[N][DIM];
double          g[N][DIM];
double          dt = .001953125;
double          Omega = .03125;
double          epsilon = .03125;
int             disallow[MAXNODES];
int             welldistant[N];
int             accepted = 0;

double
distance(double x[DIM], double y[DIM])
{
    int             j;
    double          r;
    for (r = j = 0; j < DIM; j++) {
        double          a = x[j] - y[j];
        r += a * a;
    }
    return r;
}

faraway(int i, int node)
{
    if (qdTree[node].identifier == 0) {
        return Omega * distance(x[i],
                                qdTree[node].vertex) >= qdTree[node].size;
    } else {
        return qdTree[node].identifier;
    }
}

void
qdTree_descent(int i)
{
  /**
   * As an application, we are performing a topdown search for
   * the nearest neighbors close enough to be encompassed by
   * a solid angle Omega. The main interpretation for this kind
   * of NN search is the fast calculation of the gravity forces
   * on N-body simulations.
   */
    int             level,
                    height = depth + 1;

    bzero(disallow, MAXNODES * sizeof(int));

    bzero(welldistant, N * sizeof(int));

    accepted = 0;

    for (level = 1; level < height; level++) {
        int             node;

    /** forall node in this level, do: **/
        for (node = levelbegin[level]; node < levelend[level]; node++) {
            disallow[node] = disallow[node] ||
                (faraway(i, node) && (welldistant[accepted++] = node));
        }                       // end-forall

    /** forall node in this level, do: **/
        for (node = levelbegin[level]; node < levelend[level]; node++) {
            int             q;

      /** forward node disallowance to its children **/
            for (q = 0; q < GRD; q++) {
                qdTree[node].child[q]
                    && (disallow[qdTree[node].child[q]] = disallow[node]);
            }

        }                       // end-forall

    }                           // end-for level
}

void
octree_gravity(int n)
{
    int             i,
                    k;

    for (i = 0; i < n; i++) {
        int             k;

        printf("%24.16le", mass[i]);

        for (k = 0; k < DIM; k++) {
            printf("%24.16le", x[i][k] += v[i][k] * dt);
        }

        qdTree_descent(i);

        for (bzero(g[i], sizeof(double) * DIM), k = 0; k < accepted; k++) {
            int             nu = welldistant[k],
                j;
            double          r =
                distance(qdTree[nu].vertex, x[i]) + epsilon;

            for (j = 0; j < DIM; j++) {
                double          delta = qdTree[nu].vertex[j] - x[i][j];

                g[i][j] += qdTree[nu].mass * delta / (r * sqrt(r));

            }

        }

        for (k = 0; k < DIM; k++) {
            printf("%24.16le", v[i][k] + g[i][k] * dt);
        }
#ifdef _SAVE_GRAVITY_FIELD_
        for (k = 0; k < DIM; k++) {
            printf("%24.16le", g[i][k]);
        }
#endif

        printf("\n");
    }
}

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
        fscanf(input, "%le", &mass[i]);
        Mtot += mass[i];

        for (j = 0; j < DIM; j++) {
            skipspaces(input);
            fscanf(input, "%le", &x[i][j]);
        }

        for (j = 0; j < DIM; j++) {
            skipspaces(input);
            fscanf(input, "%le", &v[i][j]);
        }

    /** seek eol **/
        while ((c = getc(input)) != '\n');

        if (skipspaces(input))
            return i + 1;
    }
}

main(int argc, char const *argv[])
{
    int             n = N,
        i;
    FILE           *input;

    if (argc > 1) {
        n = atoi(argv[1]);
        input = fopen(argv[1], "r");
    } else
        input = stdin;
    if (argc > 2) {
        Omega = atof(argv[2]);
    }
    if (argc > 3) {
        epsilon = atof(argv[3]);
    }
    if (argc > 4) {
        dt = atof(argv[4]);
    }

    qdTree_build(n = getdata(n, input));

    octree_gravity(n);

    exit(0);
}
