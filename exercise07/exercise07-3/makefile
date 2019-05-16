CC = g++
CFLAGS = -Wall -O3 --std=c++11

.PHONY: clean

exercise07-3: exercise07-3.o random.o molecular_dynamics_sim.o
	$(CC) $^ -o $@

%.o: %.cc
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm *.o exercise07-3 seed.out
