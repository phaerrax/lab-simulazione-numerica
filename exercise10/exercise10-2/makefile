CC = mpic++
CFLAGS = -Wall -O3 --std=c++11

.PHONY: clean

exercise10-2: exercise10-2.o random.o
	$(CC) $^ -o $@

%.o: %.cc
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm *.o exercise10-2 seed.out
