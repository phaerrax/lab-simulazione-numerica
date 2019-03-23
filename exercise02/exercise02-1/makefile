CC = g++
CFLAGS = -Wall -O3 --std=c++11

exercise02-1: exercise02-1.o random.o
	$(CC) random.o exercise02-1.o -o exercise02-1

exercise02-1.o: exercise02-1.cc
	$(CC) -c exercise02-1.cc -o exercise02-1.o $(CFLAGS)

random.o: random.cc random.hh
	$(CC) -c random.cc -o random.o $(CFLAGS)

clean:
	rm *.o exercise02-1 seed.out
