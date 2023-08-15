/*mksphere is an accept-reject Monte-Carlo method to generate homogeneous spherical particle distribution*/
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

unsigned short Xi[3] = { 0x1FF1, 0xAE0D, 0x3D2F };

main(int argc, char const *argv[])
{
    int n = atoi(argv[1]);
    double R = atof(argv[2]);
    double M = atof(argv[3]);
    double m = M / n;
    double u;
    double V = atof(argv[4]);
    double offset[3];
    double voffset[3];
    offset[0] = atof(argv[5]);
    offset[1] = atof(argv[6]);
    offset[2] = atof(argv[7]);
    voffset[0] = atof(argv[8]);
    voffset[1] = atof(argv[9]);
    voffset[2] = atof(argv[10]);
    u = atof(argv[11]);
    srand48(time((time_t *) Xi));
    srand48(time((time_t *) & Xi[1]));
    int i;
    for (i = 0; i < n; i++) {
        int j;
        double r, s, x[3], v[3];
        do {
            for (r = j = 0; j < 3; j++) {
                x[j] = (2 * erand48(Xi) - 1);
                r += x[j] * x[j];
            }
        } while (r > 1);
        do {
            for (s = j = 0; j < 3; j++) {
                v[j] = (2 * erand48(Xi) - 1);
                s += v[j] * v[j];
            }
        } while (s > 1);
        for (j = 0; j < 3; j++) {
            x[j] *= R;
            v[j] *= V * R * R / (r * 128 + R * R);
        }
        /*
           for (j = 0; j < 3; j++) {
           v[j] -= .25 * x[j];
           }
         */
        printf
            ("%20.16le %20.16le %20.16le %20.16le %20.16le %20.16le %20.16le %20.16le\n",
             m, x[0] + offset[0], x[1] + offset[1], x[2] + offset[2],
             v[0] + voffset[0], v[1] + voffset[1], v[2] + voffset[2], u);
    }
    return 0;
}
