#pragma once

typedef struct {
  double          d[MAXNEIGHBORINGS];
  int             list[MAXNEIGHBORINGS];
  double          radius;
  int population;
} KNN_struct;

extern KNN_struct KNN[NMAX];

typedef struct {
  double          r;
  int             left;
  int             right;
  int             ID;
} KNN_indexer;

extern int             K, KNN_population;

extern void     KNN_sort(int);

extern void     knn_search(int);

extern void     symmetric_closure(int);

