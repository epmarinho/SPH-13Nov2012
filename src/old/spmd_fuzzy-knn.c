#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <strings.h>
#include <ctype.h>
#include <spmd_octree.h>


double          mass[N];
double          x[N][DIM];

int             disallow[MAXNODES];

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

double          Omega = .03125;

isleaf(int node)
{
    int             terminal = 1,
        q;
    for (q = 0;
         q < GRD && (terminal = terminal && (qdTree[node].child[q] == 0));
         q++);
    return terminal;
}

int             fuzzyknn[N];
int             accepted = 0;

faraway(int i, int node)
{
    if (qdTree[node].population > 1) {
        double          r2 = distance(x[i], qdTree[node].vertex);
        return Omega * r2 >= qdTree[node].size;
    } else {
        fuzzyknn[accepted++] = qdTree[node].identifier;
        return 1;
    }
}

void
qdTree_descent(int i)
{
  /**
   * As an application, we are performing a topdown search for
   * the fuzzyknn neighbors close enough to be encompassed by
   * a solid angle Omega. The main interpretation for this kind
   * of NN search is the fast calculation of the gravity forces
   * on N-body simulations.
   */
    int             level;

    bzero(disallow, MAXNODES * sizeof(int));

    bzero(fuzzyknn, N * sizeof(int));

    accepted = 0;

    for (level = 1; level <= depth; level++) {
        int             node;

    /** forall node in this level, do: **/
        for (node = levelbegin[level]; node < levelend[level]; node++) {
            disallow[node] = (disallow[node] || faraway(i, node));
        }                       // end-forall

    /** forall node in this level, do: **/
        for (node = levelbegin[level]; node < levelend[level]; node++) {
            int             q;

      /** forward node disallowance to its children **/
            for (q = 0; q < GRD; q++) {
                disallow[node] && (disallow[qdTree[node].child[q]] = 1);
            }
        }                       // end-forall

    }                           // end-for level
}

void
octree_fuzzy_knn(int n)
{
    int             i = 0,
        k;

    for (i = 0; i < n; i++) {
        int             k;

        // printf("%24.16le", mass[i]);

        for (k = 0; k < DIM; k++) {
            printf("%24.16le", x[i][k]);
        }

    /** perform the topdown search for the fuzzyknn neighbors seen inside
     * the solid angle Omega */
        qdTree_descent(i);

        for (k = 0; k < accepted; k++) {
            int             nu = fuzzyknn[k],
                j;

            for (j = 0; j < DIM; j++) {
                printf("%24.16le", x[nu - 1][j]);
            }
            printf("\n");
        }
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
    for (i = 0; 1; i++) {
        int             c;

        skipspaces(input);
        fscanf(input, "%le", &mass[i]);

        for (j = 0; j < DIM; j++) {
            skipspaces(input);
            fscanf(input, "%le", &x[i][j]);
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

    qdTree_build(n = getdata(n, input));

    octree_fuzzy_knn(n);

    exit(0);
}
