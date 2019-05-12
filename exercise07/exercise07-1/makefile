CC = g++
CFLAGS = -Wall -O3 --std=c++11

.PHONY: clean

exercise07-1: exercise07-1.o metropolis_NVT.o random.o
	$(CC) $^ -o $@

warmup_analysis: warmup_analysis.o metropolis_NVT.o random.o
	$(CC) $^ -o $@

%.o: %.cc
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm *.o exercise07-1 warmup_analysis seed.out
