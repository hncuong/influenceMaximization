CC=g++
CFLAGS=-c -O3 
DEBFLAGS=-c -g 
LEADERSPATH=../leaders-3.0
OFLAGS=-DPARSEINPUT_H=1

all: SimPath

SimPath: simpath.o
	$(CC) -O3 simpath.o -o SimPath

simpath.o: simpath.cc
	$(CC) $(CFLAGS) simpath.cc	

clean:
	rm -f *.o *~