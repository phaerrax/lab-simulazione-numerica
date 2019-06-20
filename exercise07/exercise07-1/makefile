CC = g++
SRCDIR = ../..
INCLUDE = $(SRCDIR)/include
CFLAGS = -Wall -O3 --std=c++11 -I$(INCLUDE)
EXNAME = exercise07-1

.PHONY: clean

%.o: %.cc
	$(CC) $(CFLAGS) -c $< -o $@

$(EXNAME): $(EXNAME).o random.o metropolis_NVT.o
	$(CC) $^ -o $@

warmup_analysis: warmup_analysis.o metropolis_NVT.o random.o
	$(CC) $^ -o $@

%.o: %.cc
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm *.o exercise07-1 warmup_analysis seed.out
