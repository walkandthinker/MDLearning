F90=gfortran
flag=-O2 -Wall -Wextra

obj= main.o CreateAtoms.o Init.o

firstmd: ${obj}
	$(F90) $(flag) ${obj} -o firstmd

Init.o: Init.f90
	$(F90) $(flag) -c Init.f90

CreateAtoms.o: CreateAtoms.f90
	$(F90) $(flag) -c CreateAtoms.f90

main.o: main.f90
	$(F90) $(flag) -c main.f90

clean:
	rm -rf *.o *.mod firstmd