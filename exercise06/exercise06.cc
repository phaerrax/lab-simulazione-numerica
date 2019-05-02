#include <cmath>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <vector>
#include "metropolis.hh"
#include "metropolis_uniform.hh"
#include "random.hh"

std::vector<std::vector<double>> block_statistics(const std::vector<double> &, unsigned int);

int main()
{
	// Random number generator initialization
	Random rng;
	int seed[4];
	int p1, p2;
	std::ifstream Primes("Primes");
	if(Primes.is_open())
	{
		Primes >> p1 >> p2 ;
	}
	else
		std::cerr << "Unable to open Primes." << std::endl;
	Primes.close();

	std::ifstream input("seed.in");
	std::string property;
	if(input.is_open())
	{
		while(!input.eof())
		{
			input >> property;
			if(property == "RANDOMSEED")
			{
				input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
				rng.SetRandom(seed,p1,p2);
			}
		}
		input.close();
	}
	else
		std::cerr << "Unable to open seed.in." << std::endl;

	return 0;
}

std::vector<std::vector<double>> block_statistics(const std::vector<double> & x, unsigned int n_blocks)
{
	unsigned int block_size = static_cast<unsigned int>(std::round(
				static_cast<double>(x.size()) / n_blocks
				));

	// Just to make sure:
	assert(block_size * n_blocks == x.size());

	double sum(0), sum_sq(0), block_average;
	std::vector<double> row(2);
	std::vector<std::vector<double>> result;

	// - sum the values in each block;
	// - compute the average of that block;
	// - from the list of averages compute the standard dev of the mean.
	for(unsigned int i = 0; i < n_blocks; ++i)
	{
		block_average = 0;
		for(unsigned int j = 0; j < block_size; ++j)
			block_average += x[i * block_size + j];
		block_average /= block_size;
		sum           += block_average;
		sum_sq        += std::pow(block_average, 2);
		row[0]         = sum / (i + 1);
		if(i > 0)
			row[1] = std::sqrt((sum_sq / (i + 1) - std::pow(row[0], 2)) / i);
		else
			row[1] = 0;
		result.push_back(row);
	}
	return result;
}
