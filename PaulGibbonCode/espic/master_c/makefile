#
#	Makefile for 'pices' - electrostatic PIC code
#       - see file es.h for file descriptions.
#

# default flags - optimisation level 2
# - check local compiler for details ('man cc') 

#FFLAGS = -w -O2

# with debug option
FFLAGS = -w -g

# math library

LIBS = -lm

#  compiler

FC = cc

#  sources

SRCS =  pices.c init.c  loadx.c  loadv.c density.c  field.c diagnostics.c \
        plots.c  push.c  particle_bc.c histories.c add_neutral.c laser.c


#  objects

OBJS = pices.o  init.o  loadx.o loadv.o  density.o  field.o diagnostics.o \
       plots.o  push.o particle_bc.o  histories.o add_neutral.o laser.o



#  rules for compiling objects

.c.o:	
	$(FC) -c $(FFLAGS) $<

all:	pices

init.o:	init.c inputs.h

pices:   $(OBJS)
	$(FC) $(FFLAGS) $(LIBS) -o pices $(OBJS)


clean:
	rm -f basic picc *.o *~ *.data *.ps

