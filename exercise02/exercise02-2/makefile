CC = g++
CFLAGS = -Wall -O3 --std=c++11

exercise02-2: exercise02-2.o random.o
	$(CC) random.o exercise02-2.o -o exercise02-2

exercise02-2.o: exercise02-2.cc
	$(CC) -c exercise02-2.cc -o exercise02-2.o $(CFLAGS)

random.o: random.cc random.hh
	$(CC) -c random.cc -o random.o $(CFLAGS)

clean:
	rm *.o exercise02-2 seed.out
