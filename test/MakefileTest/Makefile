f90=gfortran
flag=-O3 -Wall -Wextra

OBJ=main.o PrintMod.o

test: ${OBJ}
	$(f90) $(flag) ${OBJ} -o test


PrintMod.mod: PrintMod.f90
	$(f90) $(flag) -c PrintMod.f90	
PrintMod.o: PrintMod.f90
	$(f90) $(flag) -c PrintMod.f90	


main.o:main.f90 PrintMod.mod
	$(f90) $(flag) -c main.f90

clean:
	rm *.o *.mod test
