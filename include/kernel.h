/**@<kernel.h>::*/
extern double W[NMAX][MAXNEIGHBORINGS];
extern double gradW[NMAX][MAXNEIGHBORINGS][3];
extern double h[NMAX];
extern double norm[NMAX];
extern void init_kernels(int);
extern double query2point_distance(double *, double *);
