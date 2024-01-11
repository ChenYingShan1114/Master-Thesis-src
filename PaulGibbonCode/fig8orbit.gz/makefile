FFLAGS = -w -g

singlee.o: singlee.c
	cc -c singlee.c

basicc:  singlee.o
	cc singlee.o -lm -o singec

singlefor.o: singlee.f90
	f90 -c singlee.f90 -o singlefor.o

basicf:  singlefor.o
	f90 singlefor.o -o singef

	
masterf: singlee_ma.o
	f90 singlee_ma.o -o singef

clean:
	rm singe master *.o

