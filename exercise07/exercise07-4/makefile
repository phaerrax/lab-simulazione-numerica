CC = g++
SRCDIR = ../..
INCLUDE = $(SRCDIR)/include
CFLAGS = -Wall -O3 --std=c++11 -I$(INCLUDE)
EXNAME = exercise07-4

.PHONY: clean

%.o: %.cc
	$(CC) $(CFLAGS) -c $< -o $@

$(EXNAME)_verlet: $(EXNAME)_verlet.o $(INCLUDE)/random.o molecular_dynamics_sim.o $(INCLUDE)/statistics.hh
	$(CC) $^ -o $@

$(EXNAME)_metropolis: $(EXNAME)_metropolis.o $(INCLUDE)/random.o metropolis_NVT.o $(INCLUDE)/statistics.hh
	$(CC) $^ -o $@

clean:
	rm *.o $(EXNAME) seed.out
