CC = g++
SRCDIR = ../..
INCLUDE = $(SRCDIR)/include
CFLAGS = -Wall -O3 --std=c++11 -I$(INCLUDE)
EXNAME = exercise08-2

.PHONY: clean, all
.DEFAULT_GOAL := all

%.o: %.cc
	$(CC) $(CFLAGS) -c $< -o $@

$(EXNAME): $(EXNAME).o $(INCLUDE)/random.o $(INCLUDE)/statistics.hh
	$(CC) $^ -o $@

$(EXNAME)_grid: $(EXNAME)_grid.o $(INCLUDE)/random.o $(INCLUDE)/statistics.hh
	$(CC) $^ -o $@

all: $(EXNAME) $(EXNAME)_grid

clean:
	rm *.o $(EXNAME) $(EXNAME)_grid seed.out
