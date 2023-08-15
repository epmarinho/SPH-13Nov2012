#include <strings.h>
#include <stdlib.h>
#include <math.h>
#include <constants.h>
#include <octree.h>
#include <n-body.h>

int             disallow[MAXNODES];
int             welldistant[NMAX];
double          Omega = (double) 1 / 4;
double          epsilon = (double) 1 / 64;  // it depends on the gravity
                                            // length
int             accepted;

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
    if (Octree[node].identifier == 0) {
        return Omega * distance(x[i],
                                Octree[node].vertex) >= Octree[node].size;
    } else {
        return Octree[node].identifier;
    }
}

void
Octree_descent(int i)
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

    bzero(welldistant, NMAX * sizeof(int));

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
                Octree[node].child[q]
                    && (disallow[Octree[node].child[q]] = disallow[node]);
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

        Octree_descent(i);

        for (bzero(g[i], sizeof(double) * DIM), k = 0; k < accepted; k++) {
            int             nu = welldistant[k],
                j;
            double          r = distance(Octree[nu].vertex, x[i]);

#ifdef _VARIABLE_EPSILON_
            if (Octree[nu].identifier) {
                /*
                 * now epsilon is modified by the mass influence of the
                 * gravity source in comparison to the mean particle mass
                 * of the whole n-body system
                 */
                r += epsilon * mass[Octree[nu].identifier - 1] * n / Mtot;
            }
#else
            r += epsilon;
#endif

            for (j = 0; j < DIM; j++) {
                double          delta = Octree[nu].vertex[j] - x[i][j];
                g[i][j] += Octree[nu].mass * delta / (r * sqrt(r));
            }

        }
    }
}
