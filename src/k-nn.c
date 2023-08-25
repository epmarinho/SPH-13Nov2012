#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <strings.h>
#include <ctype.h>
#include <constants.h>
#include <octree.h>
#include <macros.h>
#include <n-body.h>
#include <k-nn.h>

double          q2nd[MAXNODES];

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

KNN_struct      KNN[NMAX];

KNN_indexer     KNN_index_table[MAXNEIGHBORINGS];

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
        /*
         * visit left
         */
        if (KNN_index_table[cursor[stack]].left > 0) {
            cursor[stack + 1] = KNN_index_table[cursor[stack]].left;
            stack++;
            continue;
        }
        /*
         * visit current
         */
        if (KNN_index_table[cursor[stack]].left != -1) {
            sortlist[p++] = cursor[stack];
            KNN_index_table[cursor[stack]].left = -1;
        }
        /*
         * visit right
         */
        if (KNN_index_table[cursor[stack]].right > 0) {
            cursor[stack + 1] = KNN_index_table[cursor[stack]].right;
            KNN_index_table[cursor[stack]].right = -1;
            stack++;
            continue;
        }
        /*
         * emulate recursion return
         */
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

#ifdef _PREVIEW_SEARCH_RADIUS_
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

isfaraway(int query, int cluster)
{
    double          r;
    int             particle = Octree[cluster].identifier,
        octant;

    if (particle) {

        r = query2point_distance(x[query], x[particle - 1]);

        if (r > KNN[query].radius)
            return 1;

        if (KNN_population < K) {
            KNN_add(query, particle, r);
        } else {
            KNN_insert(query, particle, r);
        }
        return 0;

    }

    r = query2node_distance(query, cluster);

    if (r > KNN[query].radius) {
        for (octant = 0; octant < GRD; octant++) {
            int             child = Octree[cluster].child[octant];
            child && (disallow[child] = 1);
        }
        return disallow[cluster] = 1;
    }

    return 0;

}
#else
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

isfaraway(int query, int cluster)
{
    double          r;
    int             particle = Octree[cluster].identifier,
        octant;


    if (particle) {
        r = query2point_distance(x[query], x[particle - 1]);
        if (KNN_population < K) {
            KNN_add(query, particle, r);
        } else {
            if (r > KNN[query].radius)
                return 1;
            KNN_insert(query, particle, r);
        }
        return 0;
    }

    r = query2node_distance(query, cluster);
    if (r > KNN[query].radius) {
        for (octant = 0; octant < GRD; octant++) {
            int             child = Octree[cluster].child[octant];
            child && (disallow[child] = 1);
        }
        return disallow[cluster] = 1;
    }

    return 0;

}
#endif

void
knn_search(int n)
{
  /**
   * As an application, we are performing a topdown search for
   * the k-nearest neighbor according to Marinho's algorithm.
   */

    int             i,
                    level,
                    height = depth + 1,
        q,
        child,
        node;

    for (i = 0; i < n; i++) {

        bzero(q2nd, levelend[depth] * sizeof(double));
        bzero(KNN_index_table, K * sizeof(KNN_indexer));
        bzero(KNN + i, sizeof(*KNN));
        bzero(disallow, levelend[depth] * sizeof(*disallow));
        KNN_population = 0;

#ifdef _PREVIEW_SEARCH_RADIUS_
        /*
         * bottom-up find the closest ancestor with population greater than K
         */
        for (node = crossreference[i];
             node > 0 && Octree[node].population <= K;
             node = Octree[node].parent);

        /*
         * assumming 2*sigma a good choice:
         */
        KNN[i].radius = 4 * Octree[node].size;
#endif

        for (level = 1; level < height; level++) {

            for (node = levelbegin[level]; node < levelend[level]; node++) {

                disallow[node] || isfaraway(i, node);

            }

        }

        KNN[i].population = KNN_population;
    }
}

void
symmetric_closure(int n)
{
    int             i,
                    j,
                    k,
                    l;

    for (i = 0; i < n; i++) {

        for (j = 1; j < KNN[i].population; j++) {
            l = KNN[i].list[j] - 1;

            for (k = 1; k < KNN[l].population && KNN[l].list[k] != i + 1;
                 k++);

            if (k == KNN[l].population && l < MAXNEIGHBORINGS) {
                KNN[l].list[KNN[l].population++] = i + 1;
            }
        }
    }
}
