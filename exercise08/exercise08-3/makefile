CC = g++
SRCDIR = ../..
INCLUDE = $(SRCDIR)/include
CFLAGS = -Wall -O3 --std=c++11 -I$(INCLUDE)
EXNAME = exercise08-3

.PHONY: clean

%.o: %.cc
	$(CC) $(CFLAGS) -c $< -o $@

$(EXNAME): $(INCLUDE)/random.o qmc1d.o
	$(CC) $^ -o $@

clean:
	rm *.o $(EXNAME) seed.out
