CC = g++
SRCDIR = ../..
INCLUDE = $(SRCDIR)/include
CFLAGS = -Wall -O3 --std=c++11 -I$(INCLUDE)
EXNAME = exercise08-2

.PHONY: clean

%.o: %.cc
	$(CC) $(CFLAGS) -c $< -o $@

$(EXNAME): $(EXNAME).o $(INCLUDE)/random.o $(INCLUDE)/statistics.hh
	$(CC) $^ -o $@

$(EXNAME)-grid: $(EXNAME)-grid.o $(INCLUDE)/random.o $(INCLUDE)/statistics.hh
	$(CC) $^ -o $@

clean:
	rm *.o $(EXNAME) seed.out
