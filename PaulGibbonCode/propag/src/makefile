# Gnu compiler
FC = gfortran
FFLAGS =  -g -O3

# Intel compiler
#FC = ifort
#FFLAGS =  -g -fixed -O3 


.SUFFIXES: .f

SOURCES.f = absorb.f arinit.f  \
	    setup.f diagnostics.f drelp.f \
	    density.f  emhist.f \
	    emprop.f energy.f  ffts.f  \
	    flhist.f ghist.f gmovie.f \
	    grafout.f grid.f grinp.f grxy.f grzsurf.f \
	    init.f ionize.f ioutils.f launch.f \
	    outini.f powcon.f \
	    pulinit.f rezone.f spec.f \
	    stamp.f cputime.f

OBJ = $(SOURCES.f:.f=.o)

#  compile mode
.f.o:   
	$(FC) -c $(FFLAGS) $<


drelp: $(OBJ)
	$(FC) $(FFLAGS) $(OBJ) -o drelp


clean:
	rm -f drelp.o drelp core *.o


