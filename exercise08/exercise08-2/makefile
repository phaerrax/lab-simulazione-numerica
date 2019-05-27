CC = g++
CFLAGS = -Wall -O3 --std=c++11

.PHONY: clean

exercise08-2: exercise08-2.o random.o
	$(CC) $^ -o $@

exercise08-2-grid: exercise08-2-grid.o random.o
	$(CC) $^ -o $@

%.o: %.cc
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm *.o exercise08-2 exercise08-2-grid seed.out
