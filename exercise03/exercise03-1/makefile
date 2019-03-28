CC = g++
CFLAGS = -Wall -O3 --std=c++11

exercise03-1: exercise03-1.o random.o
	$(CC) random.o exercise03-1.o -o exercise03-1

exercise03-1.o: exercise03-1.cc
	$(CC) -c exercise03-1.cc -o exercise03-1.o $(CFLAGS)

random.o: random.cc random.hh
	$(CC) -c random.cc -o random.o $(CFLAGS)

clean:
	rm *.o exercise03-1 seed.out
