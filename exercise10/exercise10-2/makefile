CC = g++
MPICC = mpic++
SRCDIR = ../..
INCLUDE = $(SRCDIR)/include
CFLAGS = -Wall -O3 --std=c++11 -I$(INCLUDE)
EXNAME = exercise10-2

.PHONY: clean
.DEFAULT_GOAL := $(EXNAME)

$(EXNAME).o: $(EXNAME).cc
	$(MPICC) $(CFLAGS) -c $< -o $@

%.o: %.cc
	$(CC) $(CFLAGS) -c $< -o $@

$(EXNAME): $(EXNAME).o $(INCLUDE)/random.o
	$(MPICC) $^ -o $@

clean:
	rm *.o $(EXNAME) seed.out
