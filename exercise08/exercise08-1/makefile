CC = g++
CFLAGS = -Wall -O3 --std=c++11

.PHONY: clean

exercise08-1: exercise08-1.o random.o
	$(CC) $^ -o $@

%.o: %.cc
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm *.o exercise08-1 seed.out
