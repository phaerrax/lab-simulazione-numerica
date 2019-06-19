CC = g++
SRCDIR = ../..
INCLUDE = $(SRCDIR)/include
CFLAGS = -Wall -O3 --std=c++11 -I$(INCLUDE)

exercise01-1: exercise01-1.o random.o $(INCLUDE)/statistics.hh
	$(CC) random.o exercise01-1.o -o exercise01-1

exercise01-1.o: exercise01-1.cc $(INCLUDE)/statistics.hh
	$(CC) -c exercise01-1.cc -o exercise01-1.o $(CFLAGS)

random.o: random.cc random.hh
	$(CC) -c random.cc -o random.o $(CFLAGS)

clean:
	rm *.o exercise01-1 seed.out
