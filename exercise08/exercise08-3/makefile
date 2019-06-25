CC = g++
SRCDIR = ../..
INCLUDE = $(SRCDIR)/include
CFLAGS = -Wall -O3 --std=c++11 -I$(INCLUDE)
EXNAME = exercise08-3

.PHONY: clean, all
.DEFAULT_GOAL := all

%.o: %.cc
	$(CC) $(CFLAGS) -c $< -o $@

$(EXNAME)_bumps: $(INCLUDE)/random.o qmc1d_bumps.o
	$(CC) $^ -o $@

$(EXNAME)_flat: $(INCLUDE)/random.o qmc1d_flat.o
	$(CC) $^ -o $@

all: $(EXNAME)_bumps $(EXNAME)_flat

clean:
	rm *.o $(EXNAME)_bumps $(EXNAME)_flat seed.out
