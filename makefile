COMPILER = gfortran
CFLAGS = -O3 -Wall -g -fcheck=all -std=f2003 -pedantic 

OBJ = accuracy.o parameters.o BRENTM.o kvector.o energy.o runge_kutta.o simpson_int.o interpolation.o dgl.o matrixelements.o matrixelements_gauss.o

main: $(OBJ)
	$(COMPILER) $(CFLAGS) -o main main.f90 $(OBJ)
	make clean

BRENTM.o:
	$(COMPILER) -O3 -Wall -g -fcheck=all -pedantic -c BRENTM.f90

%.o: %.f90
	$(COMPILER) $(CFLAGS) -c $<


.PHONY: clean
.PHONY: gnuplot

clean:
	rm *.mod *.o

gnuplot:
	gnuplot -p -e "set grid;set xlabel '';set ylabel '';plot 'news.dat' using 1:3 with lines title ''"
