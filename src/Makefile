#!/bin/make -f
INCLFLAGS=-I../include
#OPTFLAGS=-floop-parallelize-all -ftree-loop-if-convert -O3
OPTFLAGS=-O3
CFLAGS=$(INCLFLAGS) $(OPTFLAGS) $(MCMODEL) $(DIRECTIVES)
LDFLAGS=-lm
#MCMODEL=-mcmodel=medium
objlist=kernel.o io-stuff.o leapfrog.o sph.o octree.o n-body.o k-nn.o
treesph: $(objlist) main.o
	$(CC) -o ../bin/$@ $^ $(LDFLAGS)
star:
	make 'DIRECTIVES=-D_INTEGRATE_ON_ENERGIES_'
indent:
	/usr/bin/indent -orig -ts2 -nuts -fca *.c
clean:
	$(RM) *.o
mostlyclean: clean
	$(RM) *~
dbg:
	make 'OPTFLAGS=-g'
gdb: dbg
set_rotation: $(objlist) set_rotation.o
	$(CC) -o ../bin/$@ $^ $(LDFLAGS)
