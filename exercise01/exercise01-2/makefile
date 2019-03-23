CC = g++
CFLAGS = -Wall -O3 --std=c++11

exercise01-2: exercise01-2.o random.o
	$(CC) random.o exercise01-2.o -o exercise01-2

exercise01-2.o: exercise01-2.cc
	$(CC) -c exercise01-2.cc -o exercise01-2.o $(CFLAGS)

random.o: random.cc random.hh
	$(CC) -c random.cc -o random.o $(CFLAGS)

clean:
	rm *.o exercise01-2 seed.out
