#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <strings.h>
#include <ctype.h>
#include <spmd_octree.h>

#define MAXNEIGHBORINGS    512
#define max(x,y)        ((x)>(y)?(x):(y))
#define min(x,y)        ((x)<(y)?(x):(y))

double          Mtot = 0;
double          mass[N];
double          x[N][DIM];
double          v[N][DIM];
double          q2nd[MAXNODES];
int             disallow[MAXNODES];

double
query2node_distance(int query, int cluster)
{
    int             k;
    double          d = 0;
    for (k = 0; k < DIM; k++) {
        double          a,
                        b;
        a = x[query][k] - Octree[Octree[cluster].parent].vertex[k];
        b = Octree[Octree[cluster].parent].vertex[k] -
            Octree[cluster].vertex[k];
        if (a * b > 0)
            d = max(d, a * a);
    }
    return q2nd[cluster] = max(d, q2nd[cluster]);
}

double
query2point_distance(double x[DIM], double y[DIM])
{
    int             j;
    double          r;
    for (r = j = 0; j < DIM; j++) {
        double          a = x[j] - y[j];
        r += a * a;
    }
    return r;
}

int             K,
                KNN_population;
struct {
    double          d[MAXNEIGHBORINGS];
    int             list[MAXNEIGHBORINGS];
    double          radius;
} KNN[N];

struct {
    double          r;
    int             left;
    int             right;
    int             ID;
} KNN_index_table[MAXNEIGHBORINGS];

extern void     KNN_sort(int);

void
KNN_add(int query, int newone, double r)
{
    int             cursor = 0;

    KNN_index_table[KNN_population].ID = newone;
    KNN[query].radius = max(KNN_index_table[KNN_population].r = r,
                            KNN[query].radius);
    if (KNN_population > 0) {
        while (1) {
            if (KNN_index_table[KNN_population].r <
                KNN_index_table[cursor].r) {
                if (KNN_index_table[cursor].left) {
                    cursor = KNN_index_table[cursor].left;
                    continue;
                }
                KNN_index_table[cursor].left = KNN_population;
            } else {
                if (KNN_index_table[cursor].right) {
                    cursor = KNN_index_table[cursor].right;
                    continue;
                }
                KNN_index_table[cursor].right = KNN_population;
            }
            break;
        }
    }
    KNN_population++;
    if (KNN_population == K) {
        KNN_sort(query);
    }
}

void
KNN_sort(int query)
{
    int             stack = 0,
        cursor[MAXNEIGHBORINGS];
    cursor[stack] = 0;
    int             sortlist[MAXNEIGHBORINGS],
                    p = 0,
        i;
    while (1) {
    /*** visit left ***/
        if (KNN_index_table[cursor[stack]].left > 0) {
            cursor[stack + 1] = KNN_index_table[cursor[stack]].left;
            stack++;
            continue;
        }
    /*** visit current ***/
        if (KNN_index_table[cursor[stack]].left != -1) {
            sortlist[p++] = cursor[stack];
            KNN_index_table[cursor[stack]].left = -1;
        }
    /*** visit right ***/
        if (KNN_index_table[cursor[stack]].right > 0) {
            cursor[stack + 1] = KNN_index_table[cursor[stack]].right;
            KNN_index_table[cursor[stack]].right = -1;
            stack++;
            continue;
        }
    /*** emulate recursion return ***/
        else if (--stack == -1) {
            break;
        } else {
            if (KNN_index_table[cursor[stack]].left > 0)
                KNN_index_table[cursor[stack]].left = 0;
        }
    }
    for (i = 0; i < K; i++) {
        int             k = sortlist[i];
        KNN[query].list[i] = KNN_index_table[k].ID;
        KNN[query].d[i] = KNN_index_table[k].r;
    }
}

void
KNN_insert(int query, int newone, double r)
{
    int             i,
                    j;
    for (i = 0; r > KNN[query].d[i]; ++i);

    for (j = KNN_population - 1; j > i; j--) {
        KNN[query].list[j] = KNN[query].list[j - 1];
        KNN[query].d[j] = KNN[query].d[j - 1];
    }
    KNN[query].list[i] = newone;
    KNN[query].d[i] = r;
    KNN[query].radius = KNN[query].d[KNN_population - 1];
}

isfaraway(int query, int cluster)
{
    double          r;
    int             particle = Octree[cluster].identifier,
        octant;

    if (particle) {

        if (KNN_population < K) {

            KNN_add(query, particle,
                    query2point_distance(x[query], x[particle - 1]));

            return 0;

        } else {

            r = query2point_distance(x[query], x[particle - 1]);

            if (r < KNN[query].radius) {

                KNN_insert(query, particle, r);

                return 0;

            } else {

                return 1;

            }
        }
    }

    r = query2node_distance(query, cluster);

    if (KNN_population < K || r <= KNN[query].radius) {

        return disallow[cluster] = 0;

    }

    for (octant = 0; octant < GRD; octant++) {
        int             child = Octree[cluster].child[octant];
        child && (disallow[child] = 1);
    }

    return disallow[cluster] = 1;

}


void
knn_search(int i, int level)
{
  /**
   * As an application, we are performing a topdown search for
   * the k-nearest neighborsaccording to Marinho's algorithm.
   */

    int
                    height = depth + 1,
        q,
        child;

    bzero(q2nd, N * sizeof(double));
    bzero(KNN_index_table, N * sizeof(*KNN_index_table));
    bzero(KNN + i, sizeof(*KNN));
    KNN_population = 0;

    for (; level < height; level++) {
        int             node;

        /*
         * forall node in this level, do: 
         */
        for (node = levelbegin[level]; node < levelend[level]; node++) {
            disallow[node] || isfaraway(i, node);
        }

#ifdef _REDUNDANT_DISALLOWANCE_PROPROPAGATION_
        /*
         * propagate node disallowance down to children: 
         */
        for (node = levelbegin[level]; node < levelend[level]; node++) {
            for (q = 0; q < GRD; q++) {
                (child = Octree[node].child[q]) &&
                    disallow[node] && (disallow[child] = 1);
            }
        }
#endif

    }
}

levellocate(int node)
{
    int             level,
                    nd;
    /*
     * locate the level the node belongs in 
     */
    for (level = 0; level < depth; level++) {
        for (nd = levelbegin[level]; nd < levelend[level]; nd++) {
            if (nd == node)
                return level;
        }
    }
    return -1;
}

void
do_knn_search(int i, int n)
{
    int               /** i, **/
                    level,
                    node,
                    nd;

  /**
  for (i = 0; i < n; i++) {
    **/
    /*
     * bottom-up find the closest ancestor with population greater than K 
     */
    for (node = crossreference[i];
         node > 0 && Octree[node].population < K;
         node = Octree[node].parent);

    /*
     * locate the level the node belongs in 
     */
    level = max(levellocate(node), 1);

    bzero(disallow, MAXNODES * sizeof(*disallow));

    for (nd = levelbegin[level]; 0 < nd && nd < levelend[level]; nd++) {
        disallow[nd] = 1;
    }

    disallow[node] = 0;

    knn_search(i, level);

    /**
  }
  **/

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
        i = 0;
    FILE           *input;

    if (argc > 1) {
        n = atoi(argv[1]);
        input = fopen(argv[1], "r");
    } else
        input = stdin;
    if (argc > 2) {
        K = atoi(argv[2]);
    }
    if (argc > 3) {
        i = atoi(argv[3]);
    }

    Octree_build(n = getdata(n, input));

    do_knn_search(i, n);

    {
        int             j,
                        k;
        for (j = 0; j < K; j++) {
            k = KNN[i].list[j] - 1;
            printf("%lg\t%lg\t%lg\t%lg\n", mass[k], x[k][0], x[k][1],
                   x[k][2]);
        }
    }

    exit(0);
}
