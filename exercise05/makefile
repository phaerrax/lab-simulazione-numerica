CC = g++
SRCDIR = ..
INCLUDE = $(SRCDIR)/include
CFLAGS = -Wall -O3 --std=c++11 -I$(INCLUDE)
EXNAME = exercise05

.PHONY: clean, all
.DEFAULT_GOAL := all

%.o: %.cc
	$(CC) $(CFLAGS) -c $< -o $@

$(EXNAME): $(EXNAME).o $(INCLUDE)/random.o $(INCLUDE)/statistics.hh
	$(CC) $^ -o $@

$(EXNAME)_raw: $(EXNAME)_raw.o $(INCLUDE)/random.o $(INCLUDE)/statistics.hh
	$(CC) $^ -o $@

all: $(EXNAME) $(EXNAME)_raw

clean:
	rm *.o $(EXNAME) $(EXNAME)_raw seed.out
