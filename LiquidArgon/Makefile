cc=g++
cflags=-O3 -Wall

objs= main.o gasdev.o initialize.o nint.o \
      verlet.o get_temp.o linearmom.o

md2:${objs}
	$(cc) $(cflags) ${objs} -o md2

nint.o:MDFunctions.h
	$(cc) $(cflags) -c nint.cpp

gasdev.o:MDFunctions.h
	$(cc) $(cflags) -c gasdev.cpp

initialize.o: MDFunctions.h
	$(cc) $(cflags) -c initialize.cpp

verlet.o: MDFunctions.h
	$(cc) $(cflags) -c verlet.cpp

linearmom.o: MDFunctions.h
	$(cc) $(cflags) -c linearmom.cpp

get_temp.o: MDFunctions.h
	$(cc) $(cflags) -c get_temp.cpp

main.o:main.cpp MDFunctions.h
	$(cc) $(cflags) -c main.cpp



clean:
	rm -rf *.o md2