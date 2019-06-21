CC = g++
SRCDIR = ../..
INCLUDE = $(SRCDIR)/include
CFLAGS = -Wall -O3 --std=c++11 -I$(INCLUDE)
EXNAME = exercise10-2

.PHONY: clean

%.o: %.cc
	$(CC) $(CFLAGS) -c $< -o $@

$(EXNAME): $(EXNAME).o $(INCLUDE)/random.o
	$(CC) $^ -o $@

clean:
	rm *.o $(EXNAME) seed.out
