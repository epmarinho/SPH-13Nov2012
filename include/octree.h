/**@<octree.h>::
 * Author: Eraldo Pereira Marinho, PhD
 * SPMD Quadtree code is a single program stream multiple data stream approach
 * to run on vector architectures like gpus etc or in clusters.
 *
 * The method consists of a level-first building of a quadtree data-structure
 * with the help of a double-hash table.
 *
 * The hash table has the form:
 * _______________________________________
 * | index | e-hash | i-hash | leaf-flag |
 * ---------------------------------------
 * |  i_1  |  eh_1  | ih_1   |    f_1    |
 * ---------------------------------------
 * |  ...  |  ...   |  ...   |   ...     |
 * |  i_n  |  eh_n  |  ih_n  |    f_n    |
 * ---------------------------------------
 * where
 * index is the entry identifier,
 * e-hash is the external hash-index which classifies entries in the same node
 * i-hash is the internal hash-index used to classify entries in the same
 * quadrant,
 * leaf-flag is an entry flag to mark the enrty for removal from the table as
 * it became a leaf.
 */
#pragma once

extern int hashTB[NMAX][4];
#define INDEX 0
#define EHASH 1
#define IHASH 2
#define LEAF  3

/*
 * Index calculation:
 *
 * i-hash = \sum_{j=1}^2{(x[i][j]-mu[j] > 0)*2^{2-j}}, where mu[] is the center
 * of mass of the node hosting the i-particle with coordinates x[].
 *
 * e-hash = 4 * e-hash + i-hash.
 */
extern int             levelbegin[MAXDEPTH],
levelend[MAXDEPTH],
depth;

extern double mean[DIM];
int ihash(int i1, int i2);
void ehash(int n);

extern struct _qdtree_ {
  /* particle identifier is a positive integer */
  int identifier;
  /* the number of particles belonging to the present quadrant */
  int population;
  /* population mass */
  double mass;
  /* node's size */
  double size;
  /* vertex is the quadrant's division point or the leaf position */
  double vertex[DIM];
  /* a cursor to the internal cluster, aka child cluster */
  int child[GRD];
  /* a parental link might be necessary */
  int parent;
} Octree[MAXNODES];

void
Octree_build(int n);


extern double          mass[NMAX];
extern double          x[NMAX][DIM];
extern double          v[NMAX][DIM];
extern double          g[NMAX][DIM];

extern  double dt;
extern int             crossreference[];
