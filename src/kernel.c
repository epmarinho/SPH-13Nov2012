/**@<kernel.c>::*/
/*
 * Smoothing kernel is modeled from the kernel based density estimation
 * theory with a compact support set by k-NN method.
 */
#include <math.h>

double
K3D_cubic_spline(double x)
{
    if (x < 1) {
        return 3.183098862e-1 * ((.75 * x - 1.5) * x * x + 1);
    }
    if (x < 2) {
        return 3.183098862e-1 * .25 * pow(2 - x, 3);
    }
    return 0;
}

double
DK3D_cubic_spline(double x)
{
    if (x < 1) {
        return 3.183098862e-1 * (2.25 * x - 3) * x;
    }
    if (x < 2) {
        return -3.183098862e-1 * .75 * pow(2 - x, 2);
    }
    return 0;
}

double
K3D_quartic_spline(double x)
{
    // x = fabs(x);
    if (x > 3)
        return 0;
    if (x > 1.5)
        return .47619047619047619047 * (.02962962962962962962 *
                                        pow(x - 3, 4));
    if (x > .5)
        return .47619047619047619047 *
            ((.3 +
              (-1.9 + (1.2 - .22962962962962962962 * x) * x) * x) * x +
             1.0875);
    return .47619047619047619047 * (1.125 +
                                    (-1 +
                                     .37037037037037037037 * x * x) * x *
                                    x);
}

double
DK3D_quartic_spline(double x)
{
    int             sgn = x < 0 ? -1 : 1;
    // x = fabs(x);
    if (x > 3)
        return 0;
    if (x > 1.5)
        return sgn * .05643738977072310403 * pow(x - 3, 3);
    if (x > .5)
        return sgn *
            (x * ((-.4373897707 * x + 1.714285714) * x - 1.809523810) +
             .1428571429);
    return x * sgn * (-.9523809524 + .7054673720 * (x * x));
}

#include <constants.h>
#include <k-nn.h>
#include <treesph.h>
#include <kernel.h>

double          W[NMAX][MAXNEIGHBORINGS];
double          gradW[NMAX][MAXNEIGHBORINGS][3];
double          h[NMAX];
double          norm[NMAX];
void
init_kernels(int N)
{
    int             i;

    knn_search(N);
#ifdef _QUARTIC_SPLINE_
    for (i = 0; i < N; h[i++] = sqrt(KNN[i].radius) / 3);
#else
    for (i = 0; i < N; h[i++] = .5 * sqrt(KNN[i].radius));
#endif
    /*
     * required for gather-scatter interpolation scheme
     */
    symmetric_closure(N);

    for (i = 0; i < N; i++) {
        int             l;
        for (l = 0; l < KNN[i].population; l++) {
            int             j = KNN[i].list[l] - 1,
                k;
            double          r_ij,
                            dx[DIM],
                            K_i,
                            K_j;

            r_ij = sqrt(query2point_distance(x[i], x[j]));
#ifdef _QUARTIC_SPLINE_
            K_i = K3D_quartic_spline(r_ij / h[i]) / pow(h[i], 3);
            K_j = K3D_quartic_spline(r_ij / h[j]) / pow(h[j], 3);
#else
            K_i = K3D_quartic_spline(r_ij / h[i]) / pow(h[i], 3);
            K_j = K3D_quartic_spline(r_ij / h[j]) / pow(h[j], 3);
#endif
            W[i][l] = .5 * (K_i + K_j);

            if (l == 0)
                continue;

            for (k = 0; k < 3; k++) {
                double          DK_i,
                                DK_j;
#ifdef _QUARTIC_SPLINE_
                DK_i = DK3D_quartic_spline(r_ij / h[i]) / pow(h[i], 4);
                DK_j = DK3D_quartic_spline(r_ij / h[j]) / pow(h[j], 4);
#else
                DK_i = DK3D_cubic_spline(r_ij / h[i]) / pow(h[i], 4);
                DK_j = DK3D_cubic_spline(r_ij / h[j]) / pow(h[j], 4);
#endif
                dx[k] = (x[i][k] - x[j][k]) / r_ij;
                gradW[i][l][k] = .5 * (DK_i + DK_j) * dx[k];
            }
        }
    }
}
