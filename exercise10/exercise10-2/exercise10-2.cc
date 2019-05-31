#include <mpi.h>
#include <iostream>
#include <functional>
#include <numeric>
#include <exception>
#include <iomanip>
#include <iterator>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <cmath>
#include <string>
#include <vector>
#include <array>
#include "random.hh"

const unsigned int n_points(30),
				   conf_elements(900);

// Define basic data structures, and the operations on them.
using point = std::array<double, 2>;
using path = std::array<point, n_points>;
point operator-(const point &);
point operator+(const point &, const point &);
point operator-(const point &, const point &);
double distance(const point &, const point &);
double path_length(const path &);
// A chromosome is a particular path, represented by an array of integers
// which indicate the order in which the cities are to be visited.
using chromosome = std::array<unsigned int, n_points>;
std::string print(const chromosome &);
// A configuration is a set of chromosomes.
using configuration = std::array<chromosome, conf_elements>;
bool is_valid(const chromosome &);

int main(int argc, char ** argv)
{ 
	MPI::Init(argc,argv);
	int size = MPI::COMM_WORLD.Get_size();
	int rank = MPI::COMM_WORLD.Get_rank();
	// Rank legend
	// 0: manager process, 
	// 1 - n-1: worker processes.
	// The manager runs the algorithm too, it's just that it has
	// alto other tasks such as the output procedures.

    // Random number generator initialization
	// ======================================
    Random rng;
    int seed[4];
    int * p = new int[size * 2];
	int p_subprocess[2];
	if(rank == 0)
	{
		std::ifstream Primes("Primes");
		std::string primes_input;
		if(Primes.is_open())
		{
			// Assign to each subprocess a different line of the file,
			// so that the generated numbers are different.
			for(int i = 0; i < size; ++i)
			{
				Primes >> p[2*i] >> p[2*i + 1];
			}
		}
		else
			std::cerr << "Unable to open Primes." << std::endl;
		Primes.close();
	}
	MPI_Scatter(
			p,              // Array of data to be sent
			2,              // Length of sent data
			MPI_INT,        // Type of sent data
			p_subprocess,   // Address to hold received data
			2,              // Length of received data
			MPI_INT,        // Type of received data
			0,              // Rank of manager process
			MPI::COMM_WORLD // Communicator
			);
	delete [] p;
	// Now each process has its own two "Primes" with that
	// it can initialise the rng with.

    std::ifstream input("seed.in");
    std::string property;
    if(input.is_open())
    {
        if(!input.eof())
        {
            input >> property;
            if(property == "RANDOMSEED")
            {
				std::cerr << "(" << rank << ") " << p_subprocess[0] << " " << p_subprocess[1] << std::endl;
                input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
                rng.SetRandom(seed, p_subprocess[0], p_subprocess[1]);
            }
        }
        input.close();
    }
    else
    {
        std::cerr << "Unable to open seed.in." << std::endl;
    }

	MPI_Barrier(MPI::COMM_WORLD);

	// Define the list of cities
	// (manager process)
	// =========================
	std::array<point, n_points> cities;
	std::string type(argv[1]);
	if(rank == 0)
	{
		if(argc != 2)
		{
			std::cerr << "Invalid input. Please specify \"circle\" or \"square\"." << std::endl;
			return 1;
		}
		unsigned int counter(0);
		while(counter < n_points)
		{
			double x, y;
			if(type == "circle")
			{
				// Generate a point on the unit circumference.
				double angle;
				angle = rng.Rannyu(0, 2 * M_PI);
				x = std::cos(angle);
				y = std::sin(angle);
			}
			else 
			{
				if(type == "square")
				{
					// Generate a point inside a square.
					x = rng.Rannyu(-1, 1);
					y = rng.Rannyu(-1, 1);
				}
				else
				{
					std::cerr << "Invalid input. Please specify \"circle\" or \"square\"." << std::endl;
					return 1;
				}
			}
			point p {x, y};
			// If there isn't already a city in that position, add it
			// to the list.
			if(std::find(cities.begin(), cities.begin() + counter, p) == cities.begin() + counter)
				cities[counter] = p;
			counter++;
		}
		std::cerr << "Manager process successfully created the list of cities" << std::endl;
	}

	// Broadcast the list of cities to other processes
	// ===============================================
	double cities_x[n_points],
		   cities_y[n_points];
	// Pack the "cities" array in two old-style arrays, one for
	// each coordinate.
	if(rank == 0)
	{
		for(unsigned int i = 0; i < n_points; ++i)
		{
			cities_x[i] = cities[i][0];
			cities_y[i] = cities[i][1];
		}
	}
	MPI_Bcast(
			cities_x,        // Pointer to the beginning of the array
			n_points,        // Length of the array
			MPI_DOUBLE,      // Type of the array values
			0,               // Manager process rank
			MPI::COMM_WORLD  // Communicator
			);
	MPI_Bcast(
			cities_y,        // Pointer to the beginning of the array
			n_points,        // Length of the array
			MPI_DOUBLE,      // Type of the array values
			0,               // Manager process rank
			MPI::COMM_WORLD  // Communicator
			);
	if(rank > 0)
	{
		// Unpack the arrays into a proper list of cities, just as the manager 
		// defined it.
		for(unsigned int i = 0; i < n_points; ++i)
		{
			cities[i][0] = cities_x[i];
			cities[i][1] = cities_y[i];
		}
		std::cerr << "(" << rank << ") city list successfully received." << std::endl;
	}
	
	// Energy function
	// ================
	// This is just the length of the path; shorter paths have less energy,
	// and the algorithm looks for the minimum of the energy.
	auto energy = [cities](const chromosome & c)
	{
		path real_path;
		for(unsigned int i = 0; i < c.size(); ++i)
			real_path[i] = cities[c[i]];
		return path_length(real_path);
	};

	// Mutations
	// =========
	// The mutations will always involve the two selected chromosomes, be they
	// replaced by their children or not.
	const unsigned int n_mutation_processes(3);
	// The mutation processes are grouped in a vector.
	std::array<std::function<chromosome (const chromosome &)>, n_mutation_processes> mutation_processes;

	// 1. Pair swap
	mutation_processes[0] = [&rng](const chromosome & c)
	{
		chromosome result(c);
		// Draw two distinct integers between 0 and n_points.
		unsigned int i1, i2;
		i1 = static_cast<unsigned int>(rng.Rannyu(0, n_points));
		do
			i2 = static_cast<unsigned int>(rng.Rannyu(0, n_points));
		while(i2 == i1);
		std::swap(result[i1], result[i2]);
		return result;
	};

	// 1. Rotation of a subsequence
	// As point 2, but rotate just of a (contiguous) subset of the chromosome.
	mutation_processes[1] = [&rng](const chromosome & c)
	{
		chromosome result(c);
		// Draw a number between 1 and n_points - 1, then rotate the array
		unsigned int seq_start = static_cast<unsigned int>(rng.Rannyu(0, n_points - 1));
		// 0 <= begin <= n_points - 2
		unsigned int seq_stop = static_cast<unsigned int>(rng.Rannyu(seq_start + 1, n_points));
		// begin + 1 <= end <= n_points - 1
		// The length of the rotated subset is always greater than 1.
		unsigned shift = static_cast<unsigned int>(rng.Rannyu(1, seq_stop - seq_start));
		std::rotate(
				result.begin() + seq_start,
				result.begin() + seq_start + shift,
				result.begin() + seq_stop
				);
		return result;
	};

	// 3. Inversion of a subsequence
	mutation_processes[2] = [&rng](const chromosome & c)
	{
		chromosome result(c);
		unsigned int begin = static_cast<unsigned int>(rng.Rannyu(0, n_points - 1));
		// 0 <= begin <= n_points - 2
		unsigned int end = static_cast<unsigned int>(rng.Rannyu(begin + 1, n_points));
		// begin + 1 <= end <= n_points - 1
		std::reverse(result.begin() + begin, result.begin() + end);
		return result;
	};

	// +----------------+
	// | Main algorithm |
	// +----------------+

	std::ofstream output_evolution("evolution_" + type + "_" + std::to_string(rank) + ".dat");

	const unsigned int col_width = 16;
	output_evolution.precision(4);
	output_evolution << std::scientific;

	chromosome current_config;

	// Generate the initial state
	// ==========================
	// Fill the chromosome with 0, 1, ..., n_points, then shuffle it.
	std::iota(current_config.begin(), current_config.end(), 0);
	// std::random_shuffle needs a function that generates (uniformly)
	// a non-negative value less than its argument.
	std::random_shuffle(
			current_config.begin(),
			current_config.end(),
			[&rng](int n) {
			return static_cast<int>(rng.Rannyu(0, n));
			}
			);

	double temperature(2 * energy(current_config));
	// Start from a temperature such that the initial acceptance rate is
	// high enough, and a lot of proposals are accepted.
	// Ex.: if it is equal to the inverse of the initial energy, then
	// the initial acceptance threshold is 1/e. A higher threshold
	// would be better: the higher the temperature, the closer the
	// ratio -energy_diff / temperature is to zero (assuming energy_diff
	// is negative, otherwise the step is always accepted).

	double cooling_rate(0.9);
	// Decrease the temperature by this factor at each step.

	unsigned int n_proposals(0),
				 accepted_proposals(0);
	double acceptance_rate(1),
		   acceptance_rate_threshold(0.025);

	unsigned int step(0);
	std::cerr << "(" << rank << ") starting algorithm." << std::endl;

	enum status_type
	{
		running = 0,
		complete = 1
	};
	status_type * status = new status_type[size];
	for(int i = 0; i < size; ++i)
		status[i] = running;
	bool all_subprocesses_complete = false;

	do
	{
		step++;

	    if(acceptance_rate != 0 && acceptance_rate < acceptance_rate_threshold)
		{
			// If the acceptance rate is already low, decreasing the temperature
			// further stalls the algorithm... I'll just terminate the program if
			// this rate drops below a certain threshold.
			// Testing for acceptance_rate == 0 avoids premature termination of
			// the loop e.g. if the very first proposal is rejected.
			
			status[rank] = complete;
			// and do nothing...
		}
		else
		{
			// Propose a new step in the Metropolis algorithm: the proposed
			// configuration is a mutation of the previous one.
			// Be careful to ensure that the detailed balance principle is
			// satisfied, e.g. the transition kernel is symmetric.
			// Anyway it shouldn't matter much since I am not really
			// trying to sample the probability function...

			// Choose at random ONE mutation and apply it.
			unsigned int index = static_cast<unsigned int>(
					rng.Rannyu(0, mutation_processes.size())
					);
			// == the index of the mutation process that was chosen.
			chromosome proposal(mutation_processes[index](current_config));
			if(!is_valid(proposal))
				throw std::runtime_error("Mutation " + std::to_string(index) + " resulted in an invalid chromosome: " + print(proposal));

			n_proposals++;

			// Test for acceptance.
			double energy_diff = energy(proposal) - energy(current_config);
			if(energy_diff < 0)
			{
				accepted_proposals++;
				current_config = proposal;
			}
			else
			{
				if(rng.Rannyu() < std::exp(-energy_diff / temperature))
				{
					accepted_proposals++;
					current_config = proposal;
				}
			}

			acceptance_rate = static_cast<double>(accepted_proposals) / n_proposals;

			if(step % 1000 == 0)
			{
				accepted_proposals = 0;
				n_proposals = 0;
				temperature *= cooling_rate;
			}
		}

		double current_energy(energy(current_config));
		double min_energy;
		// If a process exits the loop, it cannot send the data
		// anymore, and the program hangs.
		MPI_Reduce(
				&current_energy,
				&min_energy,
				1,
				MPI_DOUBLE,
				MPI_MIN,
				0,
				MPI::COMM_WORLD
				);
		if(rank == 0)
		{
			// Let the manager output the best length on a file
			output_evolution << min_energy << "\n";
		}

		// Gather the status of the subprocesses: if all of them
		// completed the algorithm, terminate the loop.
		// (status == 1 means completion.)
		status_type * status_info = new status_type[size];
		MPI_Allgather(
				&status[rank],
				1,
				MPI_INT,
				status_info,
				1,
				MPI_INT,
				MPI::COMM_WORLD
				);
		if(std::all_of(
					status_info,
					status_info + size,
					[](int i)
					{
						return i == complete;
					}
					))
		{
			// Send the signal to terminate the loop.
			all_subprocesses_complete = true;
		}
		delete [] status_info;
	}
	while(!all_subprocesses_complete);

	delete [] status;

	output_evolution.close();

	// Output procedures
	// =================
	// Output the best configuration.
	// std::ofstream best_output("best_path.dat");
	// best_output << print(current_config);
	// best_output.close();

	// Output the list of cities on a file.
	// (manager process)
	if(rank == 0)
	{
		std::ofstream cities_output_file("cities.dat");

		cities_output_file.precision(4);
		cities_output_file << std::scientific;

		for(const auto & city : cities)
		{
			for(auto x : city)
				cities_output_file << std::setw(col_width) << x;
			cities_output_file << "\n";
		}
		cities_output_file.close();
	}

	MPI::Finalize();

	return 0;
}

point operator-(const point & x)
{
	point y;
	for(unsigned int i = 0; i < x.size(); ++i)
		y[i] = -x[i];
	return y;
}

point operator+(const point & a, const point & b)
{
	point y;
	for(unsigned int i = 0; i < a.size(); ++i)
		y[i] = a[i] + b[i];
	return y;
}

point operator-(const point & a, const point & b)
{
	return a + (-b);
}

double distance(const point & a, const point & b)
{
	point diff(a - b);
	return std::sqrt(std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.));
}

double path_length(const path & p)
{
	double d(0);
	for(unsigned int i = 1; i < p.size(); ++i)
		d += distance(p[i], p[i - 1]);
	// Don't forget to add the distance from the last to the first
	// point, too.
	d += distance(p.front(), p.back());
	return d;
}

bool is_valid(const chromosome & c)
{
	// A chromosome is valid if its entries are a permutation of the
	// list (0, 1, ..., n_points).
	chromosome test;
	std::iota(test.begin(), test.end(), 0);
	return std::is_permutation(c.begin(), c.end(), test.begin());
}

std::string print(const chromosome & c)
{
	std::string s;
	for(unsigned int i = 0; i < c.size() - 1; ++i)
		s += std::to_string(c[i]) + " ";
	s += std::to_string(c.back());
	return s;
}
		
