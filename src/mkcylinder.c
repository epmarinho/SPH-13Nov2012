/*
 * mksphere is an accept-reject Monte-Carlo method to generate homogeneous 
 * spherical particle distribution
 */
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

unsigned short  Xi[3] = { 0x1FF1, 0xAE0D, 0x3D2F };

main(int argc, char const *argv[])
{
    int             n = atoi(argv[1]);
    double          R = atof(argv[2]);
    double          L = atof(argv[3]);
    double          M = atof(argv[4]);
    double          m = M / n;
    double          u;
    double          offset[3];
    double          voffset[3];
    offset[0] = atof(argv[5]);
    offset[1] = atof(argv[6]);
    offset[2] = atof(argv[7]);
    /*
     * voffset[0] = atof(argv[8]); voffset[1] = atof(argv[9]); voffset[2]
     * = atof(argv[10]); u = atof(argv[11]); 
     */
    u = atof(argv[8]);
    srand48(time((time_t *) Xi));
    srand48(time((time_t *) & Xi[1]));
    int             i;
    for (i = 0; i < n; i++) {
        int             j;
        double          r,
                        s,
                        x[3];
        do {
            for (r = j = 0; j < 2; j++) {
                x[j] = (2 * erand48(Xi) - 1);
                r += x[j] * x[j];
            }
        } while (r > 1);
        x[j] = (2 * erand48(Xi) - 1);
        for (j = 0; j < 2; j++) {
            x[j] *= R;
        }
        x[j] *= L / 2;
        printf
            ("%20.16le %20.16le %20.16le %20.16le %20.16le %20.16le %20.16le %20.16le\n",
             m, x[0] + offset[0], x[1] + offset[1], x[2] + offset[2],
             voffset[0], voffset[1], voffset[2], u);
    }
    return 0;
}
