CC = g++
CFLAGS = -Wall -O3 --std=c++11

.PHONY: clean

exercise07-2: exercise07-2.o metropolis_NVT.o random.o
	$(CC) $^ -o $@

%.o: %.cc
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm *.o exercise07-2 seed.out
