CC = g++
CFLAGS = -Wall -O3 --std=c++11

exercise01-3: exercise01-3.o random.o
	$(CC) random.o exercise01-3.o -o exercise01-3

exercise01-3.o: exercise01-3.cc
	$(CC) -c exercise01-3.cc -o exercise01-3.o $(CFLAGS)

random.o: random.cc random.hh
	$(CC) -c random.cc -o random.o $(CFLAGS)

clean:
	rm *.o exercise01-3 seed.out
