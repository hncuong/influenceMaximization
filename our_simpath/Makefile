CC=g++
CFLAGS=-c -O3 
DEBFLAGS=-c -g 
LEADERSPATH=../leaders-3.0
OFLAGS=-DPARSEINPUT_H=1

all: simpath

simpath: simpath.o
	$(CC) -O3 simpath.o -o simpath

simpath.o: simpath.cc
	$(CC) $(CFLAGS) simpath.cc simpath.hpp	

clean:
	rm -f *.o simpath *~
	rm -f temp/*.txt
