F90=gfortran
flag=-O2 -Wall -Wextra

obj= main.o CreateSample.o Init.o Force.o Integration.o

firstmd: ${obj}
	$(F90) $(flag) ${obj} -o firstmd



CreateSample.o: CreateSample.f90
	$(F90) $(flag) -c CreateSample.f90

Init.o: Init.f90
	$(F90) $(flag) -c Init.f90

Force.o: Force.f90
	$(F90) $(flag) -c Force.f90

Integration.o: Integration.f90
	$(F90) $(flag) -c Integration.f90

main.o: main.f90
	$(F90) $(flag) -c main.f90

clean:
	rm -rf *.o *.mod firstmd