#include <cmath>
#include <iomanip>
#include <algorithm>
#include <iterator>
#include <iostream>
#include <fstream>
#include <vector>
#include <cassert>
#include "random.hh"
#include "ising_1d_sim.hh"

std::vector<std::vector<double>> block_statistics(const std::vector<double> &, unsigned int);
std::vector<double> block_average(const std::vector<double> &, unsigned int);
std::vector<double> block_uncertainty(const std::vector<double> &);

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

	double in_temperature,
		   interaction_energy,
		   ext_magnetic_field;
	unsigned int n_spins,
				 metro,
				 n_blocks,
				 n_steps;

	input_parameters >> in_temperature
	                 >> n_spins
	                 >> interaction_energy
	                 >> ext_magnetic_field
	                 >> metro // if=1 Metropolis else Gibbs
	                 >> n_blocks
	                 >> n_steps;

	input_parameters.close();

	std::cout << "Temperature: "          << in_temperature     << "\n"
	          << "Number of spins: "      << n_spins            << "\n"
	          << "Exchange interaction: " << interaction_energy << "\n"
	          << "External field: "       << ext_magnetic_field << "\n";

	if(metro==1)
		std::cout << "The program perform Metropolis moves." << "\n";
	else
		std::cout << "The program perform Gibbs moves." << "\n";
	std::cout << "Number of blocks: "           << n_blocks << "\n"
	          << "Number of steps in a block: " << n_steps  << std::endl;

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

	// Define the values of the temperature that will be used to compute the
	// various physical quantities.
	std::vector<double> temperatures(50);
	double min_temperature(0.5),
		   max_temperature(2.0);
	for(unsigned int i = 0; i < temperatures.size(); ++i)
		temperatures[i] = min_temperature + (max_temperature - min_temperature) / temperatures.size() * i;

	// Define the vectors that will hold the measured values, for each value
	// of the temperature.
	std::vector<std::pair<double, double>> internal_energy,
	                                       specific_heat,
						                   magnetic_susceptibility,
										   magnetisation;

	internal_energy.reserve(temperatures.size());
	specific_heat.reserve(temperatures.size());
	magnetic_susceptibility.reserve(temperatures.size());
	magnetisation.reserve(temperatures.size());

	// Define the vectors that will hold the measured values, step by step.
	std::vector<double> energy,      // At zero external field
		                sq_energy,   // At zero external field
		                spin_sum,    // At non-zero external field
		                sq_spin_sum; // At zero external field

	// These vectors will hold the progressive values of the physical
	// quantities progressively averaged during the evolution of the
	// system, in order to perform data blocking.
	std::vector<double> avg_energy,
						avg_sq_energy,	
						avg_spin_sum,
						avg_sq_spin_sum;

	// Each time a run of the simulation is complete, from the vectors
	// defined just above we calculate the progressive (as a new block
	// of data is averaged) quantites such as the internal energy, the
	// specific heat, and so on.
	// These progressive quantities will be stored, with the respective
	// standard deviation, in the following vectors.
	std::vector<double> avg_current_internal_energy,
						std_current_internal_energy,
						avg_current_specific_heat,
						std_current_specific_heat,
						avg_current_magnetic_susceptibility,
						std_current_magnetic_susceptibility,
						avg_current_magnetisation,
						std_current_magnetisation;

	for(auto temperature : temperatures)
	{
		// Reset the vectors.
		energy.clear();
		sq_energy.clear();
		spin_sum.clear();
		sq_spin_sum.clear();
	    avg_current_internal_energy.clear();
		std_current_internal_energy.clear();
		avg_current_specific_heat.clear();
		std_current_specific_heat.clear();
		avg_current_magnetic_susceptibility.clear();
		std_current_magnetic_susceptibility.clear();
		avg_current_magnetisation.clear();
		std_current_magnetisation.clear();

		energy.reserve(n_steps);
		sq_energy.reserve(n_steps);
		spin_sum.reserve(n_steps);
		sq_spin_sum.reserve(n_steps);
	    avg_current_internal_energy.reserve(n_blocks);
		std_current_internal_energy.reserve(n_blocks);
		avg_current_specific_heat.reserve(n_blocks);
		std_current_specific_heat.reserve(n_blocks);
		avg_current_magnetic_susceptibility.reserve(n_blocks);
		std_current_magnetic_susceptibility.reserve(n_blocks);
		avg_current_magnetisation.reserve(n_blocks);
		std_current_magnetisation.reserve(n_blocks);

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
			energy.push_back(sim_zero.get_energy());
			spin_sum.push_back(sim.get_spin_sum());
			sq_spin_sum.push_back(sim_zero.get_spin_sum());

			sim.next_metropolis(rng);
			sim_zero.next_gibbs(rng);
		}

		std::transform(
				energy.begin(),
				energy.end(),
				std::back_inserter(sq_energy),
				[](double x){return x * x;}
				);

		std::transform(
				spin_sum.begin(),
				spin_sum.end(),
				std::back_inserter(sq_spin_sum),
				[](double x){return x * x;}
				);

		// Calculate average values with data blocking.
		avg_energy      = block_average(energy, n_blocks);
		avg_sq_energy   = block_average(sq_energy, n_blocks);
		avg_spin_sum    = block_average(spin_sum, n_blocks);
		avg_sq_spin_sum = block_average(sq_spin_sum, n_blocks);

		// Internal energy
		avg_current_internal_energy = avg_energy;
		std_current_internal_energy = block_uncertainty(avg_current_internal_energy);

		// Specific heat
		std::transform(
				avg_energy.begin(),
				avg_energy.end(),
				avg_sq_energy.begin(),
				std::back_inserter(avg_current_specific_heat),
				[temperature](double h, double h2){
					return (h2 - h * h) / (temperature * temperature);
				}
				);
		std_current_specific_heat = block_uncertainty(avg_current_specific_heat);

		// Magnetic susceptibility
		std::transform(
				avg_sq_spin_sum.begin(),
				avg_sq_spin_sum.end(),
				std::back_inserter(avg_current_magnetic_susceptibility),
				[temperature](double s){
					return s / temperature;
					}
				);
		std_current_magnetic_susceptibility = block_uncertainty(avg_current_magnetic_susceptibility);

		// Magnetisation
		avg_current_magnetisation = avg_spin_sum;
		std_current_magnetisation = block_uncertainty(avg_current_magnetisation);

		// Save the calculated quantities, then repeat.
		internal_energy.emplace_back(
				avg_current_internal_energy.back(),
				std_current_internal_energy.back()
				);
		specific_heat.emplace_back(
				avg_current_specific_heat.back(),
				std_current_specific_heat.back()
				);
		magnetic_susceptibility.emplace_back(
				avg_current_magnetic_susceptibility.back(),
				std_current_magnetic_susceptibility.back()
				);
		magnetisation.emplace_back(
				avg_current_magnetisation.back(),
				std_current_magnetisation.back()
				);
	}

	const unsigned int col_width = 16;

	std::ofstream output("internal_energy.dat");
	output.precision(4);
	output << std::scientific;
	for(unsigned int i = 0; i < temperatures.size(); ++i)
		output << std::setw(col_width) << temperatures[i]
		       << std::setw(col_width) << internal_energy[i].first
		       << std::setw(col_width) << internal_energy[i].second
			   << "\n";
	output.close();

	output.open("specific_heat.dat");
	for(unsigned int i = 0; i < temperatures.size(); ++i)
		output << std::setw(col_width) << temperatures[i]
		       << std::setw(col_width) << specific_heat[i].first
		       << std::setw(col_width) << specific_heat[i].second
			   << "\n";
	output.close();

	output.open("magnetic_susceptibility.dat");
	for(unsigned int i = 0; i < temperatures.size(); ++i)
		output << std::setw(col_width) << temperatures[i]
		       << std::setw(col_width) << magnetic_susceptibility[i].first
		       << std::setw(col_width) << magnetic_susceptibility[i].second
			   << "\n";
	output.close();

	output.open("magnetisation.dat");
	for(unsigned int i = 0; i < temperatures.size(); ++i)
		output << std::setw(col_width) << temperatures[i]
		       << std::setw(col_width) << magnetisation[i].first
		       << std::setw(col_width) << magnetisation[i].second
			   << "\n";
	output.close();

	return 0;
}

std::vector<double> block_average(const std::vector<double> & x, unsigned int n_blocks)
{
	// Returns a vector whose values are the averages of blocks of data
	// from x.
	std::vector<double> averages;
	averages.reserve(n_blocks);

	unsigned int block_size = static_cast<unsigned int>(std::round(
				static_cast<double>(x.size()) / n_blocks
				));

	// Just to make sure:
	assert(block_size * n_blocks == x.size());

	double block_sum;
	for(unsigned int i = 0; i < n_blocks; ++i)
	{
		block_sum = 0;
		for(unsigned int j = 0; j < block_size; ++j)
			block_sum += x[i * block_size + j];
		averages.push_back(block_sum / block_size);
	}
	return averages;
}

std::vector<double> block_uncertainty(const std::vector<double> & x)
{
	// Returns the "progressive uncertainty" of the given vector of values,
	// as if each element in the vector was a new measurement that improves
	// the uncertainty.
	// This function is best paired with block_average.
	// It's like how the uncertainty is calculated in the blocking technique
	// in Monte Carlo procedures, but here instead we already know the
	// average values.
	double avg(0), sq_avg(0);

	std::vector<double> uncertainties;
	uncertainties.reserve(x.size());

	// As usual the uncertainty on the first measurement is undefined, so
	// we start from the second one.
	uncertainties.push_back(0);
	unsigned int count(2);
	
	for(auto v : x)
	{
		avg    += v;
		sq_avg += v * v;
		uncertainties.push_back(
				std::sqrt((sq_avg / count - std::pow(avg / count, 2)) / (count - 1))
				);
		++count;
	}

	return uncertainties;
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
