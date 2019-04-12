CC = g++
CFLAGS = -Wall -O3 --std=c++11

.PHONY: clean

exercise05: exercise05.o random.o
	$(CC) $^ -o $@

%.o: %.cc
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm *.o exercise05 seed.out
