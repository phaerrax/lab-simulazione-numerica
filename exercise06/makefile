CC = g++
CFLAGS = -Wall -O3 --std=c++11

.PHONY: clean

exercise06: exercise06.o random.o
	$(CC) $^ -o $@

%.o: %.cc
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm *.o exercise06 seed.out
