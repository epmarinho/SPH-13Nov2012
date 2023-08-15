#include <rasterfile.h>

#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#include <libgen.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>

void
setcolor(double x, char *r, char *g, char *b ){
    if (0<=x && x<0.125){
        r[0]=0xDC;
        g[0]=0xDC;
        b[0]=0xDC;
    }else if (0.125<=x && x<0.250){
        r[0]=0xFF;
        g[0]=0x07;
        b[0]=0x01;
    }else if (0.250<=x && x<0.475){
        r[0]=0xFF;
        g[0]=0x65;
        b[0]=0x01;
    }else if (0.475<=x && x<0.500){
        r[0]=0xFF;
        g[0]=0xFF;
        b[0]=0x00;
    }else if (0.500<=x && x<0.625){
        r[0]=0xFF;
        g[0]=0xFF;
        b[0]=0xFF;
    }else if (0.750<=x && x<0.875){
        r[0]=0x00;
        g[0]=0x20;
        b[0]=0xFE;
    }else if (0.875<=x && x<1){
        r[0]=0x66;
        g[0]=0x00;
        b[0]=0xFF;
    }else if (x<=1){
        r[0]=0x7F;
        g[0]=0x01;
        b[0]=0x7F;
    }
}
/*
 * Get device definitions
 */
/*
 * Device definitions
 */
/*
 * Raster constants
 */
#define COLORDEPTH 24

typedef struct rasterfile RASTER;

char           *pixmap;

RASTER         *header;

/*
 * End of Device definitions
 */

/** image parameters **/
char           *pixmap;

int             NX /** image width in pixels **/ ,
NY /** image height in pixels **/ ,
NXB /** image width in bytes **/ ;

/** data frame: x in [a:b], y in [c:d] **/
double          a,
b,
c,
d;

char *R,*G,*B;

#define max(x,y) (x>y?x:y)
#define min(x,y) (x<y?x:y)

/*
 *
 * Programa de conversao de palavra entre maquinas tipo SUN e DEC-Alpha.
 *
 * A sequencia de bytes numa SUN Sparc XX e' inversa a das DECs
 *
 */

invert_word(int word)
{
    char            byte,
    *p;

    register int    k;

    int             i;

    int             l = sizeof(int);  /* recebe o comprimento da palavra */

    /*
     * Opcionalmente:
     */
    /*
     * int l = sizeof(word);
     */

    p = (char *) &word;           /* toma o endereco de word na pilha */

    for (k = 0; k < l - k; k++) {
        byte = *(p + k);
        *(p + k) = *(p + l - k - 1);
        *(p + l - k - 1) = byte;
    }

    return *((int *) p);
}

void
SaveRaster(int fd)
{
    write(fd, (void *) header, sizeof(struct rasterfile));
    write(fd, (void *) pixmap, NXB * NY);
    close(fd);
}

void
SetupRasHeader(void)
{
    int             i,
    *byte;
    /**  allocate memory space for the header structure **/
    header = (RASTER *) calloc(sizeof(struct rasterfile), 1);
    /** setup header parameters **/
    header->ras_magic = RAS_MAGIC;  /* Magic number */
    header->ras_width = NX;       /* horizontal size */
    header->ras_height = NY;      /* vertical size */
    header->ras_depth = COLORDEPTH; /* color depth in bits */
    header->ras_length = (NXB = COLORDEPTH / 8 * NX + (NX % 2)) * NY; /* total
    * number
    * of
    * bytes */
    header->ras_type = RT_STANDARD;
    header->ras_maptype = RMT_NONE;
    header->ras_maplength = 0;
    /** Allocate pixmap space in bytes **/
    pixmap = calloc(header->ras_length, 1);
    #ifdef LITTLE_ENDIAN
    byte = (int *) header;
    for (i = 0; i < sizeof(struct rasterfile); i += 4)
        *byte++ = invert_word(*byte);
    #endif
}

double         *x,
*y, *u;
int            *i,
*j;
int             xcol = 1,
ycol = 2,vxcol = 4,vycol = 5;

void
plotXY(int n)
{
    int             l;
    for (l = 0; l < n; l++) {
        setcolor(u[l],R+l,G+l,B+l);
        /** BGR in little-endian archs, and RGB in big-endian **/

        /** blue  **/pixmap[3 * j[l] + 0 + i[l] * NXB] = B[l];
        /** green **/pixmap[3 * j[l] + 1 + i[l] * NXB] = G[l];
        /** red   **/pixmap[3 * j[l] + 2 + i[l] * NXB] = R[l];
    }
}

/** data reading **/
read_data(int n, FILE * input_data)
{
    int             l,k;
    double vx, vy, vmax=0, vmin=1e38;
    for (l = 0; !feof(input_data); l++) {
        int             c,
        j;

        if (xcol < ycol) {

            for (j = 0; j < xcol; j++) {
                while (isspace(c = getc(input_data)));
                ungetc(c, input_data);
                fscanf(input_data, "%lg", x + l);
            }

            for (; j < ycol; j++) {
                while (isspace(c = getc(input_data)));
                ungetc(c, input_data);
                fscanf(input_data, "%lg", y + l);
            }

            if(vxcol < vycol){

                for (; j < vxcol; j++) {
                    while (isspace(c = getc(input_data)));
                    ungetc(c, input_data);
                    fscanf(input_data, "%lg", &vx);
                    u[l]=vx*vx;
                }

                for (; j < vycol; j++) {
                    while (isspace(c = getc(input_data)));
                    ungetc(c, input_data);
                    fscanf(input_data, "%lg", &vy);
                    u[l]+=vy*vy;
                    vmax = max(vmax,u[l]);
                    vmin = min(vmin,u[l]);
                }

            } else {

                for (; j < vycol; j++) {
                    while (isspace(c = getc(input_data)));
                    ungetc(c, input_data);
                    fscanf(input_data, "%lg", &vy);
                    u[l]+=vy*vy;
                }

                for (; j < vxcol; j++) {
                    while (isspace(c = getc(input_data)));
                    ungetc(c, input_data);
                    fscanf(input_data, "%lg", &vx);
                    u[l]+=vx*vx;
                    vmax = max(vmax,u[l]);
                    vmin = min(vmin,u[l]);
                }

            }

        } else {

            for (j = 0; j < ycol; j++) {
                while (isspace(c = getc(input_data)));
                ungetc(c, input_data);
                fscanf(input_data, "%lg", y + l);
            }

            for (; j < xcol; j++) {
                while (isspace(c = getc(input_data)));
                ungetc(c, input_data);
                fscanf(input_data, "%lg", x + l);
            }

            if(vxcol < vycol){

                for (; j < vxcol; j++) {
                    while (isspace(c = getc(input_data)));
                    ungetc(c, input_data);
                    fscanf(input_data, "%lg", &vx);
                    u[l]=vx*vx;
                }

                for (; j < vycol; j++) {
                    while (isspace(c = getc(input_data)));
                    ungetc(c, input_data);
                    fscanf(input_data, "%lg", &vy);
                    u[l]+=vy*vy;
                    vmax = max(vmax,u[l]);
                    vmin = min(vmin,u[l]);
                }

            } else {

                for (; j < vycol; j++) {
                    while (isspace(c = getc(input_data)));
                    ungetc(c, input_data);
                    fscanf(input_data, "%lg", &vy);
                    u[l]+=vy*vy;
                }

                for (; j < vxcol; j++) {
                    while (isspace(c = getc(input_data)));
                    ungetc(c, input_data);
                    fscanf(input_data, "%lg", &vx);
                    u[l]+=vx*vx;
                    vmax = max(vmax,u[l]);
                    vmin = min(vmin,u[l]);
                }

            }
        }
        while (!((c = getc(input_data)) == '\n' || c == -1));
    }

    for(k = 0; k < l; k++){
        u[k] = (u[k] - vmin) / (vmax - vmin);
    }

    return l;
}

XY2JI(int n)
{
    int             l,
    k = 0;
    double          sigma = (d - c) / (b - a);
    NY = sigma * NX;
    for (l = 0; l < n; l++) {
        int             jj,
        ii;
        jj = NX * (x[l] - a) / (b - a);
        ii = NY * (1 - (y[l] - c) / (d - c));
        if (0 <= jj && jj < NX && 0 <= ii && ii < NY) {
            j[k] = jj;
            i[k] = ii;
            k++;
        }
    }
    return k;
}

/** main part **/
main(int argc, char const *argv[])
{
    if (argc < 11) {
        fprintf(stderr, "11 parameters required... exiting\n");
        exit(-2);
    } else {
        int             parm = 1,
        n = 0;

        FILE           *buffer = fopen(argv[parm++], "r");

        if (buffer == NULL) {
            fprintf(stderr, "cannot read %s\n", argv[1]);
            exit(-1);
        } else {
            int             c;
            while ((c = getc(buffer)) != -1) {
                if (c == '\n')
                    n++;
            }
            rewind(buffer);
            x = malloc(n * sizeof(double));
            y = malloc(n * sizeof(double));
            u = malloc(n * sizeof(double));
            R = malloc(n * sizeof(double));
            G = malloc(n * sizeof(double));
            B = malloc(n * sizeof(double));
            j = malloc(n * sizeof(double));
            i = malloc(n * sizeof(double));
        }

        a = atof(argv[parm++]);
        b = atof(argv[parm++]);
        c = atof(argv[parm++]);
        d = atof(argv[parm++]);
        xcol = atoi(argv[parm++]);
        ycol = atoi(argv[parm++]);
        vxcol = atoi(argv[parm++]);
        vycol = atoi(argv[parm++]);
        NX = atoi(argv[parm++]);

        read_data(n, buffer);

        n = XY2JI(n);

        SetupRasHeader();

        plotXY(n);

        SaveRaster(1);

        exit(0);

    }
}
