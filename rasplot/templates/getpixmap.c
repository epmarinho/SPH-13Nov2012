/** This is a sample program to extract pixel map from a SUN Rasterfile **/
/** Author: Eraldo Pereira Marinho
    Changes:  2008, May 30
              2008, Jun 02
 **/

#define max(x,y) (x>y?x:y)
#define min(x,y) (x<y?x:y)

#include "rasterfile.h"

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

invertlong(int x)
{
  int             y;
  int             j;
  for (j = 0; j < 4; j++) {
    ((char *) &y)[3 - j] = ((char *) &x)[j];
  }
  return y;
}


main(int argc, char const *argv[])
{
  int             NX,
                  NY,
                  Npixels;
  size_t          linesiz;
  char            cmdname[32];
  struct rasterfile rasheader;
  int             nr;
  unsigned char **r,
                **g,
                **b;
  unsigned      **grey;
  unsigned        mingrey = 1 << 31,
      maxgrey = 0,
      minred = 1 << 31,
      maxred = 0,
      mingreen = 1 << 31,
      maxgreen = 0,
      minblue = 1 << 31,
      maxblue = 0;
  char           *room;
  double          sigr,
                  sigg,
                  sigb,
                  siggrey;
  double          mred = 0,
      mgreen = 0,
      mblue = 0,
      mgrey,
      msqred = 0,
      msqgreen = 0,
      msqblue = 0,
      msqgrey = 0;
  strcpy(cmdname, basename((char *) argv[0]));
  /*
   * acquire the rasterfile path and allocate a file descriptor for it
   * -- otherwise, assume stdin 
   */
  if (argc > 1) {
    int             fd = open(argv[1], O_RDONLY);
    if (fd < 0) {
      char           *errmsg = ": no such rasterfile: ";
      int             errcode = errno;
      write(2, cmdname, strlen(cmdname));
      write(2, errmsg, strlen(errmsg));
      write(2, argv[1], strlen(argv[1]));
      write(2, "\n", 1);
      exit(errcode);
    }
    dup2(fd, 0);
  }
  if ((nr = read(0, &rasheader, sizeof(rasheader))) == 0) {
    char           *errmsg = ": rasterfile is empty\n";
    write(2, cmdname, strlen(cmdname));
    write(2, errmsg, strlen(errmsg));
    exit(-1);
  }
  /*
   * get the virtual image width in pixels 
   */
  NX = invertlong(rasheader.ras_width);
  /*
   * get the image height 
   */
  NY = invertlong(rasheader.ras_height);
  /*
   * total number of pixels
   */
  Npixels = NX * NY;
  /*
   * allocate grid space
   */
  r = (unsigned char **) malloc(NY * sizeof(void));
  g = (unsigned char **) malloc(NY * sizeof(void));
  b = (unsigned char **) malloc(NY * sizeof(void));
  grey = (unsigned **) malloc(NY * sizeof(void));
  /*
   * get the image line size in bytes (multiple of 16 bits) 
   */
  linesiz = invertlong(rasheader.ras_length) / NY;
  /*
   * allocate room space for image line 
   */
  room = malloc(linesiz);
  /*
   * fill up a RGB grid 
   */
  {
    int             i;
    for (i = 0; i < NY; i++) {
      /*
       * read an entire image line 
       */
      int             nr = read(0, room, linesiz),
          j,
          k;
      if (nr < linesiz) {
        fprintf(stderr, "%s: rasterfile is corrupt\n", cmdname);
        exit(-2);
      }
      r[i] = (unsigned char *) malloc(NX);
      g[i] = (unsigned char *) malloc(NX);
      b[i] = (unsigned char *) malloc(NX);
      grey[i] = (unsigned *) malloc(NX * sizeof(int));
      for (k = j = 0; k < NX; (j += 3), (k++)) {
        /*
         * image_raw_data 
         */
        b[i][k] = room[j + 0];
        g[i][k] = room[j + 1];
        r[i][k] = room[j + 2];
        grey[i][k] = b[i][k] + g[i][k] + r[i][k];
        /*
         * image_statistics 
         */
        /** extreme values **/
        mingrey = min(mingrey, grey[i][k]);
        minred = min(minred, r[i][k]);
        mingreen = min(mingreen, g[i][k]);
        minblue = min(minblue, b[i][k]);
        maxgrey = max(maxgrey, grey[i][k]);
        maxred = max(maxred, r[i][k]);
        maxgreen = max(maxgreen, g[i][k]);
        maxblue = max(maxblue, b[i][k]);
        /** mean values **/
        /*** first momentum */
        mred += r[i][k];
        mgreen += g[i][k];
        mblue += b[i][k];
        mgrey += grey[i][k];
        /*** second momentum ***/
        msqred += r[i][k] * r[i][k];
        msqgreen += g[i][k] * g[i][k];
        msqblue += b[i][k] * b[i][k];
        msqgrey += grey[i][k] * grey[i][k];
      }
    }
    /*
     * mean RGB+grey colors 
     */
    mred /= Npixels;
    mgreen /= Npixels;
    mblue /= Npixels;
    mgrey /= Npixels;
    msqred /= Npixels;
    msqgreen /= Npixels;
    msqblue /= Npixels;
    msqgrey /= Npixels;
    /*
     * color standard deviation 
     */
    sigr = sqrt(msqred - mred * mred);
    sigg = sqrt(msqgreen - mgreen * mgreen);
    sigb = sqrt(msqblue - mblue * mblue);
    siggrey = sqrt(msqgrey - mgrey * mgrey);
  }
  printf("%s: pixel map was succefully extracted from rasterfile ",
         cmdname);
  if (argc > 1)
    printf("%s", argv[1]);
  printf("\n");
  return 0;
}
