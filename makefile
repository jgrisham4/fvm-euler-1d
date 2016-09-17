FC=gfortran
FCFLAGS=-O2 
SRC_FILES = $(wildcard *.f90)

TARGETS = test_mesh test_solver

all : $(TARGETS)

%.o : %.f90
	$(FC) -c $(FCFLAGS) $<

test_mesh : mesh.o test_mesh.o
	$(FC) $(FCFLAGS) -o $@ $^

test_solver : mesh.o utils.o limiters.o flux.o solver.o test_solver.o
	$(FC) $(FCFLAGS) -o $@ $^

clean :
	rm -rf *.o *.tec *.dat *.out *.mod $(TARGETS)
