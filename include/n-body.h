extern int             disallow[MAXNODES];
extern int             welldistant[NMAX];
extern double          Omega;
extern double          epsilon;
extern int             accepted;

extern double Mtot;

extern double distance(double x[DIM], double y[DIM]);
extern faraway(int i, int node);
extern void Octree_descent(int i);
extern void octree_gravity(int n);