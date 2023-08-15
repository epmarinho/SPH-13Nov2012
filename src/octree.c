#include <strings.h>

#include <constants.h>
#include <octree.h>

int             levelbegin[MAXDEPTH],
levelend[MAXDEPTH],
depth;

int             node = 0,
next_child_node = 1;

/* The quadtree data structure is modeled here: */
struct _qdtree_ Octree[MAXNODES];

/* It is sometimes interesting to have a cross-reference index to locate the node
 * which a given particle is stored */
int             crossreference[NMAX];

/* Author: Eraldo Pereira Marinho, PhD
 * Reference: Marinho 1997, Doctoral Thesis, Chapter 8.
 *
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
 * ---------------------------------------
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

int             hashTB[NMAX][4];
/*
 * Index calculation:
 *
 * i-hash = \sum_{j=1}^2{(x[i][j]-mu[j] > 0)*2^{2-j}}, where mu[] is the center
 * of mass of the node hosting the i-particle with coordinates x[].
 *
 * e-hash = 4 * e-hash + i-hash.
 */

int             collision[GRD];

ihash(int i1, int i2)
{
  /* the [i1, i2)-range is semantically a node in the current
   * quadtree level. */

  if (i2 == i1 + 1) {
    int             j,
    k = hashTB[i1][INDEX] - 1;
    hashTB[i1][EHASH] = hashTB[i1][IHASH] = -1;
    hashTB[i1][LEAF] = 1;
    Octree[crossreference[k] = node].identifier = hashTB[i1][INDEX];
    Octree[node].mass = mass[k];
    for (j = 0; j < DIM; j++) {
      Octree[node].vertex[j] = x[k][j];
    }
    Octree[node].size = Octree[Octree[node].parent].size;
    return 1;

  } else {

    int             i;
    double          Mtot = 0;

    bzero(Octree[node].vertex, sizeof(double) * DIM);

    bzero(collision, sizeof(int) * GRD);

    Octree[node].size = 0;

    /* for parallel */
    for (i = i1; i < i2; i++) {
      int             j,
      k = hashTB[i][INDEX] - 1;
      for (j = 0; j < DIM; j++) {
        Octree[node].vertex[j] += mass[k] * x[k][j];
        Octree[node].size += mass[k] * x[k][j] * x[k][j];
      }
      Mtot += mass[k];
    }

    Octree[node].size /= Mtot;

    /* for parallel */
    {
      int             j;
      for (j = 0; j < DIM; j++) {
        Octree[node].vertex[j] /= Mtot;
        Octree[node].size -=
        Octree[node].vertex[j] * Octree[node].vertex[j];
      }
    }

    Octree[node].mass = Mtot;

    /* for parallel */
    for (i = i1; i < i2; i++) {
      int             j,
      k = hashTB[i][INDEX] - 1;
      for (hashTB[i][IHASH] = j = 0; j < DIM; j++) {
        hashTB[i][IHASH] *= 2;
        hashTB[i][IHASH] += (x[k][j] - Octree[node].vertex[j] > 0);
      }
    }

    for (i = i1; i < i2; i++) {
      collision[hashTB[i][IHASH]]++;
    }

    /*
     * duplicates check
     */
    {
      int             has_singularity;
      for (has_singularity = i = 0; i < GRD; i++) {
        has_singularity += (collision[i] == 0);
      }
      has_singularity = (has_singularity == GRD - 1);
    }

    bzero(Octree[node].child, sizeof(int) * GRD);

    for (i = 0; i < GRD; i++) {
      if (collision[i] > 0) {
        Octree[next_child_node].parent = node;
        Octree[node].child[i] = next_child_node++;
      }
    }

    return i2 - i1;
  }
}

void
ehash(int n)
{
  int             i,
  newindex;
  for (i = 0; i < n; i++) {
    hashTB[i][EHASH] *= GRD;
    hashTB[i][EHASH] += hashTB[i][IHASH];
  }

  for (i = 0; i < n - 1; i++) {
    int             j;
    for (j = i + 1; j < n; j++) {
      if (hashTB[i][EHASH] > hashTB[j][EHASH]) {
        int             k;
        for (k = 0; k < 4; k++) {
          int             aux = hashTB[i][k];
          hashTB[i][k] = hashTB[j][k];
          hashTB[j][k] = aux;
        }
      }
    }
  }

  for (newindex = 0, i = 1; i < n; i++) {
    if (hashTB[i - 1][EHASH] < hashTB[i][EHASH])
      hashTB[i - 1][EHASH] = newindex++;
    else
      hashTB[i - 1][EHASH] = newindex;
  }
  hashTB[i - 1][EHASH] = newindex;

}

leafprune(int n)
{
  int             i,
  j = 0;
  int             TBaux[NMAX][4];
  for (i = 0; i < n; i++) {
    if (hashTB[i][LEAF] == 0) {
      int             k;
      for (k = 0; k < 4; k++) {
        TBaux[j][k] = hashTB[i][k];
      }
      j++;
    }
  }
  if (j == n)
    return n;
  for (i = 0; i < j; i++) {
    int             k;
    for (k = 0; k < 4; k++) {
      hashTB[i][k] = TBaux[i][k];
    }
  }
  return j;
}

void
hashTBinit(int n)
{
  int             i;
  for (i = 0; i < n; i++)
    hashTB[i][INDEX] = i + 1;
}

#include <stdlib.h>
#include <stdio.h>

#define max(x,y)        (x>y?x:y)
#define min(x,y)        (x<y?x:y)

void
Octree_build(int n)
{
  int             level;

  hashTBinit(n);

  for (level = 0; n; level++) {

    int             i1,
    i2;

    levelbegin[level] = node;
    levelend[level] = next_child_node;

    /* foreach node, do */
    for (i1 = 0; i1 < n; node++, i1 = i2) {

      for (i2 = i1 + 1;
           i2 < n && hashTB[i2 - 1][EHASH] == hashTB[i2][EHASH];
      i2++);

      Octree[node].population = ihash(i1, i2);

    }
    /* end-foreach */

    /* delete all just formed leafs */
    n = leafprune(n);

    if (n == 0)
      break;

    ehash(n);
  }

  depth = level;
}
