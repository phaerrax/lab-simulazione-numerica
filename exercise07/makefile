CC = g++
CFLAGS = -Wall -O3 --std=c++11

.PHONY: clean

exercise07: exercise07.o metropolis_NVT.o random.o
	$(CC) $^ -o $@

%.o: %.cc
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm *.o exercise07 seed.out
