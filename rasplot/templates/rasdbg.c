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

invert32(int x)
{
  int             y;
  int             j;
  for (j = 0; j < 4; j++) {
    ((char *) &y)[3 - j] = ((char *) &x)[j];
  }
  return y;
}

/*
 * This is a simple code to display header content from Sun Rasterfiles,
 * which is quite useful on developing low-level image processing tools 
 */
main(int argc, char const *argv[])
{
  char            cmdname[32];
  struct rasterfile rasheader;
  int             nr;
  strcpy(cmdname, argv[0]);
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
  nr = read(0, &rasheader, sizeof(rasheader));
  if (nr == 0) {
    char           *errmsg = ": rasterfile is empty\n";
    write(2, cmdname, strlen(cmdname));
    write(2, errmsg, strlen(errmsg));
    exit(-1);
  }
  /*
   * reordering for 68000 endian convention 
   */
  rasheader.ras_magic = invert32(rasheader.ras_magic);
  rasheader.ras_width = invert32(rasheader.ras_width);
  rasheader.ras_height = invert32(rasheader.ras_height);
  rasheader.ras_depth = invert32(rasheader.ras_depth);
  rasheader.ras_length = invert32(rasheader.ras_length);
  rasheader.ras_type = invert32(rasheader.ras_type);
  rasheader.ras_maptype = invert32(rasheader.ras_maptype);
  rasheader.ras_maplength = invert32(rasheader.ras_maplength);

  printf("rasterfile header:\n");
  printf("magic number (should be 0x59a66a95): 0x%8x\n", rasheader.ras_magic);  /* magic 
                                                                                 * number 
                                                                                 */
  printf("width: %i\n", rasheader.ras_width);   /* width (pixels) of image 
                                                 */
  printf("height: %i\n", rasheader.ras_height); /* height (pixels) of
                                                 * image */
  printf("depth: %i\n", rasheader.ras_depth);   /* depth (1, 8, or 24
                                                 * bits) of pixel */
  printf("length: %i\n", rasheader.ras_length); /* length (bytes) of image 
                                                 */
  printf("type: 0x%x\n", rasheader.ras_type);   /* type of file; see RT_*
                                                 * below */
  printf("maptype: 0x%x\n", rasheader.ras_maptype);     /* type of
                                                         * colormap; see
                                                         * RMT_* below */
  printf("maplength: 0x%x\n", rasheader.ras_maplength); /* length (bytes)
                                                         * * * of
                                                         * following * * * 
                                                         * * * * * map */
  printf("checksum:\n");
  printf("line room (bytes): %i\n",
         rasheader.ras_length / rasheader.ras_height);
  printf("line size (bytes): %i\n", 3 * rasheader.ras_width);
  return 0;
}
