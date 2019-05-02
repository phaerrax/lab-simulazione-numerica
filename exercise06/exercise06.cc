#include <cmath>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <vector>
#include <cassert>
#include "random.hh"
#include "ising_1d_sim.hh"

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

	// Initialisation procedure
	// ========================
	std::cout << "Monte Carlo simulation of a classic 1D Ising model,\n"
	          << "using nearest neighbour interaction.\n"
	          << "The program uses k_B=1 and mu_B=1." << std::endl;
	// Read the input parameters from a file.
	// They have to be given in a very specific way...
	std::string input_parameters_file("input.dat");
	std::cout << "Read input parameters from " << input_parameters_file << "." << std::endl;
	std::ifstream input_parameters(input_parameters_file);

	double temperature,
		   interaction_energy,
		   ext_magnetic_field;
	unsigned int n_spins,
				 metro,
				 n_blocks,
				 n_steps;

	input_parameters >> temperature
	                 >> n_spins
	                 >> interaction_energy
	                 >> ext_magnetic_field
	                 >> metro // if=1 Metropolis else Gibbs
	                 >> n_blocks
	                 >> n_steps;

	input_parameters.close();

	std::cout << "Temperature: "          << temperature        << "\n"
	          << "Number of spins: "      << n_spins            << "\n"
	          << "Exchange interaction: " << interaction_energy << "\n"
	          << "External field: "       << ext_magnetic_field << "\n";

	if(metro==1)
		std::cout << "The program perform Metropolis moves." << "\n";
	else
		std::cout << "The program perform Gibbs moves." << "\n";
	std::cout << "Number of blocks: "           << n_blocks << "\n"
	          << "Number of steps in a block: " << n_steps  << std::endl;

	// Prepare arrays for measurements.
	unsigned int n_props(4), // Number of observables
	             iu(0),      // Energy
	             ic(1),      // Heat capacity
	             im(2),      // Magnetization
	             ix(3);      // Magnetic susceptibility

	// Initialize the vector that will contain the measured physical quantities
	// at each step.
	std::vector<double> istantaneous_quantities(4);

	// Once the system has reached equilibrium, we can start taking averages.
	// At each new step, the system evolves to a new state, according to the
	// Metropolis or the Gibbs algorithm.
	// The physical quantities of interest are
	// - internal energy,
	// - specific heat,
	// - magnetic susceptibility,
	// - magnetisation;
	// in order to calculate them we will need the average value of the
	// Hamiltonian, of its square, and the square of the average of the sum
	// of all spins, at zero external magnetic field, and the average of
	// the sum of all spins at a positive magnetic field.
	// Therefore, at each step, we need to know the Hamiltonian and the
	// sum of all spins for both zero and positive magnetic field.
	//
	// Since the Hamiltonian depends on the external magnetic field, we have
	// to consider two different system at a time, one with the field set to
	// zero and the other to the desired positive value (but initialised in
	// the same way), since the evolution in the space of configuration is
	// different whether the field is zero or not. 
	//
	// The uncertainties for these four quantities will be estimated using a
	// blocking technique: each time a new step is reached, the quantity is
	// calculated again; eventually we will have, for each quantity, a list
	// of values, from which we can extract the average and the standard
	// deviation.

	// Generate an initial configuration.
	ising_1d_sim sim(n_spins, rng);

	// Set initial parameters.
	sim.set_inverse_temperature(1. / temperature);
	sim.set_interaction_energy(interaction_energy);
	sim.set_external_magnetic_field(ext_magnetic_field);

	// Define a new model with the external field set to zero.
	// We will use sim to calculate the magnetisation, sim_zero to calculate
	// the other quantities.
	ising_1d_sim sim_zero(sim);
	sim_zero.set_external_magnetic_field(0);

	// Start the simulation.
	for(unsigned int step = 0; step < n_steps; ++step)
	{
		sim.next_metropolis(rng);
		sim_zero.next_gibbs(rng);
	}
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
