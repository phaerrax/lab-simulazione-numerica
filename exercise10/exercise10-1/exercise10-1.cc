#include <iostream>
#include <functional>
#include <numeric>
#include <exception>
#include <iomanip>
#include <iterator>
#include <algorithm>
#include <fstream>
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
    {
        std::cerr << "Unable to open seed.in." << std::endl;
    }

	// Define the list of cities
	// =========================
    if(argc != 2)
    {
        std::cerr << "Invalid input. Please specify \"circle\" or \"square\"." << std::endl;
        return 1;
    }
    std::string type(argv[1]);
	std::vector<point> cities;
	unsigned int generated_cities(0);
	while(generated_cities <= 30)
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
		if(std::find(cities.begin(), cities.end(), p) == cities.end())
			cities.push_back(p);
		generated_cities++;
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

	unsigned int n_proposals(0),
				 accepted_proposals(0);

	std::ofstream output_evolution(type + "/evolution.dat");
	
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

	double acceptance_rate,
		   acceptance_rate_threshold(0.01);

	unsigned int step(0);
	do
	{
		step++;
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

			// Output the current configuration, and its energy (= its length).
			output_evolution << energy(current_config)
				             << " " << print(current_config) << "\n";
			std::cerr << "\rTemperature: " << temperature
				      << " / acceptance rate: " << acceptance_rate;
		}
					
		// If the acceptance rate is already low, decreasing the temperature
		// further stalls the algorithm... I'll just terminate the program if
		// this rate drops below a certain threshold.
		// Testing for current_temp_acceptance_rate == 0 avoids premature
		// termination of the loop e.g. if the very first proposal is rejected.
	}
	while(acceptance_rate == 0 || acceptance_rate > acceptance_rate_threshold);

	std::cerr << std::endl;

	output_evolution.close();

	// Output procedures
	// =================
	// Output the best configuration.
	std::ofstream best_output(type + "/best_path.dat");
	best_output << print(current_config);
	best_output.close();

	// Output the list of cities on a file.
	std::ofstream cities_output_file(type + "/cities.dat");

    cities_output_file.precision(4);
    cities_output_file << std::scientific;

	for(const auto & city : cities)
	{
		for(auto x : city)
			cities_output_file << std::setw(col_width) << x;
		cities_output_file << "\n";
	}
	cities_output_file.close();

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
		
