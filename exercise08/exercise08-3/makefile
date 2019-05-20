CC = g++
CFLAGS = -Wall -O3 --std=c++11

.PHONY: clean

qmc1d: qmc1d.o random.o
	$(CC) $^ -o $@

%.o: %.cc
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm *.o qmc1d seed.out
