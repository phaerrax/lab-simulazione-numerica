CC = g++
CFLAGS = -Wall -O3 --std=c++11

.PHONY: clean

exercise09: exercise09.o random.o
	$(CC) $^ -o $@

%.o: %.cc
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm *.o exercise09 seed.out
