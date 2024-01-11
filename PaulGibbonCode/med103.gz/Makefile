OBJ = med103.o 
FFLAGS = -w -O2 -fixed
#  debug mode
#FFLAGS = -g
FC = ifort

med103: med103.o
	$(FC) -o med103 $(OBJ)
med103.o: med103.f
	$(FC) $(FFLAGS) -c med103.f
debug:  debug.o
	$(FC) $(OBJ) -o debug
debug.o:  med103.f
	$(FC)  -C -w -c -g med103.f
medgraf.o:  medgraf.f 
ioutils.o:  ioutils.f
	$(FC) $(FFLAGS) -c ioutils.f

medgraf:	medgraf.o ioutils.o
	$(FC) medgraf.o ioutils.o -o medgraf

timehist.o:  timehist.f 
timehist:	timehist.o 
	$(FC) timehist.o ioutils.o -o timehist

ionhist.o:  ionhist.f 
ionhist:	ionhist.o 
	$(FC) ionhist.o ioutils.o -o ionhist

default: all
all: 	med103 medgraf timehist ionhist

clean:
	rm -f med103.o med103 debug debug.o medgraf.o medgraf

