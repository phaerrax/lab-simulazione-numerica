CC = g++
SRCDIR = ../..
INCLUDE = $(SRCDIR)/include
CFLAGS = -Wall -O3 --std=c++11 -I$(INCLUDE)
EXNAME = exercise07-4

.PHONY: clean, all
.DEFAULT_GOAL := all

%.o: %.cc
	$(CC) $(CFLAGS) -c $< -o $@

$(EXNAME)_verlet: $(EXNAME)_verlet.o $(INCLUDE)/random.o molecular_dynamics_sim.o $(INCLUDE)/statistics.hh
	$(CC) $^ -o $@

$(EXNAME)_metropolis: $(EXNAME)_metropolis.o $(INCLUDE)/random.o metropolis_NVT.o $(INCLUDE)/statistics.hh
	$(CC) $^ -o $@

all: $(EXNAME)_verlet $(EXNAME)_metropolis

clean:
	rm *.o $(EXNAME)_verlet $(EXNAME)_metropolis seed.out
