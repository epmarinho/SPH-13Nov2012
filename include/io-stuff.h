#pragma once
extern double          Mtot;
extern double          mass[NMAX];
extern double          x[NMAX][DIM];
extern double          v[NMAX][DIM];
extern double          g[NMAX][DIM];
extern double          dt;

extern skipspaces(FILE * );
extern getdata(int n, FILE * );
extern void show_result(int n, FILE *output);
