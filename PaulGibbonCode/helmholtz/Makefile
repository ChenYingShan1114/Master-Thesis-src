# Makefile for Helmholtz solver
F90 = ifort 
.SUFFIXES: .f90
FFLAGS = -g

.f90.o:
	$(F90) -c $(FFLAGS) $<

helm:	helmholtz.o
	$(F90) -o helm helmholtz.o 

brunel: brunel.o
	$(F90) -o brunel brunel.o 

clean:
	/bin/rm -f *.o helm
	echo ..cleaned 
