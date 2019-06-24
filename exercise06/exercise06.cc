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

/*
   These two functions are the "unpackaged" versions of the usual
   block_statistics function I used in other exercises. Since in
   the program there will be some quantities for which the
   average value and the uncertainty cannot be calculated
   simultaneously (at least, not without distorting the
   block_statistics function), I need the separate functions.
*/
template <class ForwardIterator, class OutputIterator>
void block_average(ForwardIterator, ForwardIterator, OutputIterator, unsigned int);
template <class ForwardIterator, class OutputIterator>
void block_uncertainty(ForwardIterator, ForwardIterator, OutputIterator);

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

	double equilibration_steps,
		   interaction_energy,
		   ext_magnetic_field;
	unsigned int n_spins,
				 metro,
				 block_size,
				 n_steps;

	input_parameters >> equilibration_steps
	                 >> n_spins
	                 >> interaction_energy
	                 >> ext_magnetic_field
	                 >> metro // if == 1 Metropolis else Gibbs
	                 >> block_size
	                 >> n_steps;

	input_parameters.close();

	std::cout << "Equilibration steps: "  << equilibration_steps << "\n"
	          << "Number of spins: "      << n_spins             << "\n"
	          << "Exchange interaction: " << interaction_energy  << "\n"
	          << "External field: "       << ext_magnetic_field  << "\n";

	if(metro == 1)
		std::cout << "The program perform Metropolis moves." << "\n";
	else
		std::cout << "The program perform Gibbs moves."       << "\n";

	std::cout << "Block size: "      << block_size << "\n"
	          << "Number of steps: " << n_steps    << std::endl;

	/*
	   Once the system has reached equilibrium, we can start taking averages.
	   At each new step, the system evolves to a new state, according to the
	   Metropolis or the Gibbs algorithm.
	   The physical quantities of interest are
	   - internal energy,
	   - specific heat,
	   - magnetic susceptibility,
	   - magnetisation;
	   in order to calculate them we will need the average value of the
	   Hamiltonian, of its square, and the square of the average of the sum
	   of all spins, at zero external magnetic field, and the average of
	   the sum of all spins at a positive magnetic field.
	   Therefore, at each step, we need to know the Hamiltonian and the
	   sum of all spins for both zero and positive magnetic field.

	   Since the Hamiltonian depends on the external magnetic field, we have
	   to consider two different system at a time, one with the field set to
	   zero and the other to the desired positive value (but initialised in
	   the same way), since the evolution in the space of configuration is
	   different whether the field is zero or not.

	   The uncertainties for these four quantities will be estimated using a
	   blocking technique: each time a new step is reached, the quantity is
	   calculated again; eventually we will have, for each quantity, a list
	   of values, from which we can extract the average and the standard
	   deviation.
	*/

	// Define the values of the temperature that will be used to compute the
	// various physical quantities.
	std::vector<double> temperatures(50);
	double min_temperature(0.5),
		   max_temperature(2.0);
	for(unsigned int i = 0; i < temperatures.size(); ++i)
		temperatures[i] = min_temperature + (max_temperature - min_temperature) / temperatures.size() * i;

	/*
	   Define the vectors that will hold the measured values, for each value
	   of the temperature.
	   Each sub-vector represents a different quantity:
	   [0] internal energy;
	   [1] specific heat;
	   [2] magnetic susceptibility;
	   [3] magnetisation.
	*/
	std::vector<std::vector<std::pair<double, double>>> measurements(4);
	std::vector<std::string> measurement_names = {
		"internal_energy",
		"specific_heat",
		"magnetic_susceptibility",
		"magnetisation"
	};

	for(auto & m : measurements)
		m.reserve(temperatures.size());

	// Define the vectors that will hold the measured values, step by step.
	std::vector<double> energy,
						// At zero external field, needed for the internal
						// energy and the specific heat.
		                sq_energy,
						// At zero external field, needed for the
						// specific heat.
		                spin_sum,
						// At non-zero external field, needed for the
						// magnetisation.
		                spin_sum_zero,
						// At zero external field, needed for the
						// magnetic susceptibility.
		                sq_spin_sum;
						// At zero external field, needed for the
						// magnetic susceptility (will hold the square
						// of the values in spin_sum_zero).

	// These vectors will hold the progressive of the physical
	// quantities, progressively averaged during the evolution of the
	// system, in order to perform data blocking.
	std::vector<double> avg_energy,
						avg_sq_energy,	
						avg_spin_sum,
						avg_sq_spin_sum;

	/*
	   Each time a run of the simulation is complete, from the vectors
	   defined just above we calculate the progressive (as a new block
	   of data is averaged) quantites such as the internal energy, the
	   specific heat, and so on.
	   These progressive quantities will be stored, with the respective
	   standard deviation, in the following vectors.
	*/
	std::vector<std::vector<double>> avg_current_measurements(measurements.size()),
						             std_current_measurements(measurements.size());

	for(auto temperature : temperatures)
	{
		// Reset the vectors.
		energy.clear();
		sq_energy.clear();
		spin_sum.clear();
		spin_sum_zero.clear();
		sq_spin_sum.clear();
		energy.reserve(n_steps);
		sq_energy.reserve(n_steps);
		spin_sum.reserve(n_steps);
		spin_sum_zero.reserve(n_steps);
		sq_spin_sum.reserve(n_steps);

	    for(auto & m : avg_current_measurements)
			m.clear();
		for(auto & m : std_current_measurements)
			m.clear();

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
		for(unsigned int step = 0; step < equilibration_steps; ++step)
		{
			if(metro == 1)
			{
				sim.next_metropolis(rng);
				sim_zero.next_metropolis(rng);
			}
			else
			{
				sim.next_gibbs(rng);
				sim_zero.next_gibbs(rng);
			}
		}
		for(unsigned int step = equilibration_steps; step < n_steps; ++step)
		{
			energy.push_back(sim_zero.get_energy());
			spin_sum.push_back(sim.get_spin_sum());
			spin_sum_zero.push_back(sim_zero.get_spin_sum());

			if(metro == 1)
			{
				sim.next_metropolis(rng);
				sim_zero.next_metropolis(rng);
			}
			else
			{
				sim.next_gibbs(rng);
				sim_zero.next_gibbs(rng);
			}
		}

		// Calculate the square of the Hamiltonian and the sum of spins.
		std::transform(
				std::begin(energy),
				std::end(energy),
				std::back_inserter(sq_energy),
				[](double x){return x * x;}
				);

		std::transform(
				std::begin(spin_sum_zero),
				std::end(spin_sum_zero),
				std::back_inserter(sq_spin_sum),
				// The square of the sum of all spins is needed only
				// in the h = 0 case, so I will not indicate it
				// in the variable name anymore.
				[](double x){return x * x;}
				);

		// Calculate average values of the Hamiltonian and the sum of spins
		// using data blocking... 
		block_average(
				std::begin(energy),
				std::end(energy),
				std::back_inserter(avg_energy),
				block_size
				);
		block_average(
				std::begin(sq_energy),
				std::end(sq_energy),
				std::back_inserter(avg_sq_energy),
				block_size
				);
		block_average(
				std::begin(spin_sum),
				std::end(spin_sum),
				std::back_inserter(avg_spin_sum),
				block_size
				);
		block_average(
				std::begin(sq_spin_sum),
				std::end(sq_spin_sum),
				std::back_inserter(avg_sq_spin_sum),
				block_size
				);

		// ...and from these values, the progressive average values of the
		// physical quantities of interest:
		// [0] internal energy
		avg_current_measurements[0] = avg_energy;

		// [1] specific heat
		std::transform(
				avg_energy.begin(),
				avg_energy.end(),
				avg_sq_energy.begin(),
				std::back_inserter(avg_current_measurements[1]),
				[temperature](double h, double h2){
					return (h2 - std::pow(h, 2)) * std::pow(temperature, -2);
				}
				);

		// [2] magnetic susceptibility
		std::transform(
				avg_sq_spin_sum.begin(),
				avg_sq_spin_sum.end(),
				std::back_inserter(avg_current_measurements[2]),
				[temperature](double s){
					return s / temperature;
					}
				);

		// [3] magnetisation
		avg_current_measurements[3] = avg_spin_sum;

		// From the list of average values, compute the progressive
		// uncertainties.
		for(unsigned int i = 0; i < avg_current_measurements.size(); ++i)
			block_uncertainty(
					std::begin(avg_current_measurements[i]),
					std::end(avg_current_measurements[i]),
					std::back_inserter(std_current_measurements[i])
					);

		// Save the calculated quantities, then repeat.
		// (The last value of each measurement is taken as the "final" value
		// of the relative quantity at the current temperature.)
		for(unsigned int i = 0; i < avg_current_measurements.size(); ++i)
			measurements[i].emplace_back(
					avg_current_measurements[i].back(),
					std_current_measurements[i].back()
					);
	}

	const unsigned int col_width = 16;

	std::ofstream output;
	output.precision(4);
	output << std::scientific;

	for(unsigned int i = 0; i < measurements.size(); ++i)
	{
		std::string filename(measurement_names[i]);
		if(metro == 1)
			filename += "_metropolis.dat";
		else
			filename += "_gibbs.dat";
		output.open(filename);
		for(unsigned int j = 0; j < temperatures.size(); ++j)
			output << std::setw(col_width) << temperatures[j]
				   << std::setw(col_width) << measurements[i][j].first
				   << std::setw(col_width) << measurements[i][j].second
				   << "\n";
		output.close();
	}

	rng.SaveSeed();
	return 0;
}

template <class ForwardIterator, class OutputIterator>
void block_average(ForwardIterator first, ForwardIterator last, OutputIterator avg_first, unsigned int block_size)
{
    // If there are not enough values to form a block, typically
    // when we are near the end of the input data, the remaining
    // the points are discarded.
    double sum(0), block_average(0);
    for(unsigned int blocks = 1; std::distance(first, last) >= block_size; ++blocks)
    {
        // Sum the values in each block.
        block_average = 0;
        for(unsigned int j = 0; j < block_size; ++j)
        {
            block_average += *first;
            ++first;
            // At the end of each iteration of the main loop the first
            // iterator will have advanced by block_size positions.
        }

        // Compute the average of that block.
        block_average /= block_size;
        sum           += block_average;

        // Save the results.
        *avg_first     = sum / blocks;
        ++avg_first;
        // If the output iterators are back_inserters or similar iterators,
        // they advance ad the next position when they are assigned a value,
        // and the operator++ has no effect on them.
    }

	return;
}

template <class ForwardIterator, class OutputIterator>
void block_uncertainty(ForwardIterator first, ForwardIterator last, OutputIterator err_first)
{
	// Returns the "progressive uncertainty" of the given vector of values,
	// as if each element in the vector was a new measurement that improves
	// the uncertainty.
	// This function is best paired with block_average.
	// It's like how the uncertainty is calculated in the blocking technique
	// in Monte Carlo procedures, but here instead we already know the
	// average values.
	double avg(0), sq_avg(0);

	// As usual the uncertainty on the first measurement is undefined, so
	// we start from the second one.
	*err_first = 0;
	err_first++;

	unsigned int count(2);
	
	while(first != last)
	{
		avg    += *first;
		sq_avg += std::pow(*first, 2);
		*err_first = std::sqrt((sq_avg / count - std::pow(avg / count, 2)) / (count - 1));
		++count;
		++err_first;
		++first;
	}

	return;
}
