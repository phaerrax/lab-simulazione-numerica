CC = g++
CFLAGS = -Wall -O3 --std=c++11

.PHONY: clean

exercise10-1: exercise10-1.o random.o
	$(CC) $^ -o $@

%.o: %.cc
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm *.o exercise10-1 seed.out
