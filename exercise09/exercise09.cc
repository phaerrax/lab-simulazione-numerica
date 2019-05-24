#include <iostream>
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
// A configuration is a set of chromosomes.
using configuration = std::array<chromosome, conf_elements>;
bool is_valid(const chromosome &);
std::string to_string(const chromosome &);

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

	// Define the list of cities
	// =========================
	std::vector<point> cities;
	unsigned int generated_cities(0);
	while(generated_cities <= 30)
	{
		// Generate a point on the unit circumference.
		double x, y, angle;
		angle = rng.Rannyu(0, 2 * M_PI);
		x = std::cos(angle);
		y = std::sin(angle);
		point p {x, y};
		// If there isn't already a city in that position, add it
		// to the list.
		if(std::find(cities.begin(), cities.end(), p) == cities.end())
			cities.push_back(p);
		generated_cities++;
	}
	
	// Fitness function
	// ================
	// The higher the value of the function (to be
	// evaluated on a chromosome) the better the fitness.
	auto fitness = [cities](const chromosome & c)
	{
		path real_path;
		for(unsigned int i = 0; i < c.size(); ++i)
			real_path[i] = cities[c[i]];
		// A chromosome is fitter the shorter is the path it represents.
		return -path_length(real_path);
	};

	// Mutations
	// =========
	// The mutations will always involve the two selected chromosomes, be they
	// replaced by their children or not.
	const unsigned int n_mutation_processes(4);
	// The mutation processes are grouped in a vector.
	std::array<std::function<void (chromosome &)>, n_mutation_processes> mutation_processes;
	// This vector contains the probability of each mutation process to occur.
	std::array<double, n_mutation_processes> mutation_probability;
	for(auto & p : mutation_probability)
		p = 0.01;

	// 1. Pair swap
	mutation_processes[0] = [&rng](chromosome & c)
	{
		// Draw two distinct integers between 0 and n_points.
		unsigned int i1, i2;
		i1 = static_cast<unsigned int>(rng.Rannyu(0, n_points));
		do
			i2 = static_cast<unsigned int>(rng.Rannyu(0, n_points));
		while(i2 == i1);
		std::swap(c[i1], c[i2]);
		return;
	};

	// 2. Rotation
	mutation_processes[1] = [&rng](chromosome & c)
	{
		// Draw a number between 1 and n_points - 1, then rotate the array
		// that number of times.
		unsigned int shift = static_cast<unsigned int>(rng.Rannyu(1, n_points));
		std::rotate(
				c.begin(),
				c.begin() + shift,
				c.end()
				);
		return;
	};

	// 3. Rotation of a subsequence
	// As point 2, but rotate just of a (contiguous) subset of the chromosome.
	mutation_processes[2] = [&rng](chromosome & c)
	{
		// Draw a number between 1 and n_points - 1, then rotate the array
		unsigned int seq_start = static_cast<unsigned int>(rng.Rannyu(0, n_points - 1));
		// 0 <= begin <= n_points - 2
		unsigned int seq_stop = static_cast<unsigned int>(rng.Rannyu(seq_start + 1, n_points));
		// begin + 1 <= end <= n_points - 1
		// The length of the rotated subset is always greater than 1.
		unsigned shift = static_cast<unsigned int>(rng.Rannyu(1, seq_stop - seq_start));
		std::rotate(
				c.begin() + seq_start,
				c.begin() + seq_start + shift,
				c.begin() + seq_stop
				);
		return;
	};

	// 4. Inversion of a subsequence
	mutation_processes[3] = [&rng](chromosome & c)
	{
		unsigned int begin = static_cast<unsigned int>(rng.Rannyu(0, n_points - 1));
		// 0 <= begin <= n_points - 2
		unsigned int end = static_cast<unsigned int>(rng.Rannyu(begin + 1, n_points));
		// begin + 1 <= end <= n_points - 1
		std::reverse(c.begin() + begin, c.begin() + end);
		return;
	};

	// +----------------+
	// | Main algorithm |
	// +----------------+
	
	// Generate the initial population
	// ===============================
	std::array<chromosome, 900> chromosomes;
	for(auto & x : chromosomes)
	{
		// Fill the chromosome with 0, 1, ..., n_points, then shuffle it.
		std::iota(x.begin(), x.end(), 0);
		// std::random_shuffle needs a function that generates (uniformly)
		// a non-negative value less than its argument.
		std::random_shuffle(
				x.begin(),
				x.end(),
				[&rng](int n) {
					return static_cast<int>(rng.Rannyu(0, n));
					}
				);
	}

	// Sort initially the configuration in order of increasing fitness.
	std::sort(
			chromosomes.begin(),
			chromosomes.end(),
			[fitness](const chromosome & a, const chromosome & b) {
				return (fitness(a) < fitness(b));
			}
			);

	chromosome best_fit, previous_best_fit(chromosomes.back());
	// Hold the shortest path at each step, and its previous step.
	double crossover_probability(0.8);
	unsigned int steps_without_improvements(0),
	             max_steps_without_improvements(100),
				 generation(0);
	// If the algorithm is "stuck" for more than the given number 
	// of steps, then exit and return the best fit.
	// steps_without_improvements is incremented if the best_fit has not
	// changed from the previous step.
	
	std::vector<chromosome> new_generation;
	// Holds the newly generated chromosomes, at each iteration of the loop,
	// before they are copied in the configuration array.

	// Define a function that draws an integer from 0 to conf_elements, such that
	// higher values are favoured.
	auto skewed_draw = [&rng]()
	// (No need to capture conf_elements, since it's not an automatic variable.
	// It compiles fine on g++, but it may give an error on other compilers.)
	{
		// If x is uniformly distributed in [0, 1), then the pdf of x^p (p > 0) is
		// 1/p * x^(1/p - 1), therefore higher values are more likely to be drawn
		// if I choose p < 1. The more p is close to 1, the more the skewed
		// distribution is close to the uniform one; the more p is close to 0,
		// the more the distribution is peaked at 1.
		return static_cast<unsigned int>(conf_elements * std::pow(rng.Rannyu(), 0.5));
	};

	do
	{
		generation++;
		std::cerr << "\rGeneration #" << generation;

		new_generation.clear();

		while(new_generation.size() < conf_elements)
		{
			// Parents selection
			// =================
			// Shorter paths are in the last positions of the list of configurations,
			// therefore the selection algorithm should favour high numbers, so that
			// fitter chromosomes are more eligible to the "reproduction" process.

			// Draw two (distinct) chromosomes from the list.
			unsigned int parent1, parent2;
			parent1 = skewed_draw();
			do
				parent2 = skewed_draw();
			while(parent1 == parent2);

			chromosome child1, child2;

			// Crossover
			// =========
			// Roll for crossover.
			if(rng.Rannyu() < crossover_probability)
			{
				// Select (randomly) a point at which to cut the chromosomes.
				unsigned int cut_point = static_cast<unsigned int>(rng.Rannyu(0, n_points));
				// Define two new chromosomes, and copy the parents' information
				// up to the cut point.
				std::copy(
						chromosomes[parent1].begin(),
						chromosomes[parent1].begin() + cut_point,
						child1.begin()
						);
				std::copy(
						chromosomes[parent2].begin(),
						chromosomes[parent2].begin() + cut_point,
						child2.begin()
						);
				// Complete the chromosomes with the missing cities, adding them in
				// the order in which they appear in the "other parent" (child1 with
				// parent2, child2 with parent1).
				std::copy_if(
						chromosomes[parent1].begin(),
						chromosomes[parent1].end(),
						child2.begin() + cut_point,
						[&child2, cut_point](unsigned int n)
						{
							return std::find(child2.begin(), child2.begin() + cut_point, n) == child2.begin() + cut_point;
						}
						);
				std::copy_if(
						chromosomes[parent2].begin(),
						chromosomes[parent2].end(),
						child1.begin() + cut_point,
						[&child1, cut_point](unsigned int n)
						{
							return std::find(child1.begin(), child1.begin() + cut_point, n) == child1.begin() + cut_point;
						}
						);
			}
			else
			{
				child1 = chromosomes[parent1];
				child2 = chromosomes[parent2];
			}

			// Mutate the two newly generated chromosomes (or their parents, if
			// crossover didn't happen).
			for(unsigned int i = 0; i < mutation_processes.size(); ++i)
			{
				if(rng.Rannyu() < mutation_probability[i])
					mutation_processes[i](child1);
				if(!is_valid(child1))
					throw std::runtime_error("Mutation " + std::to_string(i) + " resulted in an invalid chromosome: " + to_string(child1));
				if(rng.Rannyu() < mutation_probability[i])
					mutation_processes[i](child2);
				if(!is_valid(child2))
					throw std::runtime_error("Mutation " + std::to_string(i) + " resulted in an invalid chromosome: " + to_string(child1));
			}

			// The two new chromosomes are ready.
			new_generation.push_back(std::move(child1));
			new_generation.push_back(std::move(child2));
		}

		// The new generation is now complete, so we can replace the old
		// chromosomes with it.
		std::copy(
				new_generation.begin(),
				new_generation.end(),
				chromosomes.begin()
				);

		// Sort again the configuration in order of increasing fitness.
		// Now, the last chromosome is the candidate for the minimum
		// distance path.
		std::sort(
				chromosomes.begin(),
				chromosomes.end(),
				[fitness](const chromosome & a, const chromosome & b)
				{
					return fitness(a) < fitness(b);
				}
				);

		previous_best_fit = best_fit;
		best_fit = chromosomes.back();

		if(best_fit == previous_best_fit)
			steps_without_improvements++;
	}
	while(steps_without_improvements < max_steps_without_improvements);

	// Output procedures
	// =================
	// Output the list of cities on a file.
	std::ofstream cities_output_file("cities.dat");

    const unsigned int col_width = 16;
    cities_output_file.precision(4);
    cities_output_file << std::scientific;

	for(const auto & city : cities)
	{
		for(auto x : city)
			cities_output_file << std::setw(col_width) << x;
		cities_output_file << "\n";
	}
	cities_output_file.close();

	// Output the minimum configuration on a file.
	std::ofstream minimum_output_file("minimum.dat");
	std::ostream_iterator<unsigned int> minimum_output(minimum_output_file, "\n");
	std::copy(best_fit.begin(), best_fit.end(), minimum_output);
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

std::string to_string(const chromosome & c)
{
	std::string s;
	for(unsigned int i = 0; i < c.size() - 1; ++i)
		s += std::to_string(c[i]) + " ";
	s += std::to_string(c.back());
	return s;
}
		
