#include <cmath>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <iterator>
#include <numeric>
#include <functional>
#include <vector>
#include "metropolis.hh"
#include "metropolis_uniform.hh"
#include "random.hh"

std::vector<std::vector<double>> block_statistics(const std::vector<double> &, unsigned int);

int main(int argc, char * argv[])
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

	if(argc != 6)
	{
		std::cerr << "Invalid input.\n Syntax: " << argv[0] << " mu_min mu_max sigma_min sigma_max #points." << std::endl;
		return 1;
	}

    const unsigned int dim = 1;

	// Command line arguments
	double mu_min(std::stod(argv[1])),
	       mu_max(std::stod(argv[2])),
	       sigma_min(std::stod(argv[3])),
	       sigma_max(std::stod(argv[4]));

	unsigned int n_grid_points(std::stoi(argv[5]));

	std::vector<double> mu_values(n_grid_points),
						sigma_values(n_grid_points);
	for(unsigned int i = 0; i < mu_values.size(); ++i)
		mu_values[i] = mu_min + (mu_max - mu_min) / mu_values.size() * i;
	for(unsigned int i = 0; i < sigma_values.size(); ++i)
		sigma_values[i] = sigma_min + (sigma_max - sigma_min) / sigma_values.size() * i;

    double L(1.2);
    metropolis_uniform<dim> metro(-L, L);
    
    // Output the progressive values of the average energy and its standard
    // deviation, obtained with a blocking technique, to a file.
	std::ofstream output_file("grid.dat");
    const unsigned int col_width = 16;
    output_file.precision(4);
    output_file << std::scientific;

	unsigned int n_blocks(100);
    std::array<double, dim> start {0};
	unsigned int calc(0); // To show some form of progress.

	for(auto mu : mu_values)
	{
		for(auto sigma : sigma_values)
		{
			metro.set_starting_point(start);

			auto trial_wf = [mu, sigma](std::array<double, dim> r) {
				double x = r[0];
				return std::exp(-0.5 * std::pow((x - mu) / sigma, 2)) + std::exp(-0.5 * std::pow((x + mu) / sigma, 2));
			};

			auto pdf = [mu, sigma, trial_wf](std::array<double, dim> r) {
				return std::pow(trial_wf(r), 2);
			};

			auto trial_wf_2nd = [mu, sigma](std::array<double, dim> r) {
				double x = r[0];
				double plus  = std::exp(-0.5 * std::pow((x - mu) / sigma, 2));
				double minus = std::exp(-0.5 * std::pow((x + mu) / sigma, 2));
				return -0.5 * std::pow(sigma, -2) * (
						(plus + minus) +
						std::pow(sigma, -2) * (std::pow(x - mu, 2) * minus + std::pow(x + mu, 2) * plus)
						);
			};

			auto potential_energy = [](std::array<double, dim> r) {
				double x = r[0];
				return std::pow(x, 4) - 5. / 2. * std::pow(x, 2);
			};

			auto energy = [mu, sigma, trial_wf, trial_wf_2nd, potential_energy](std::array<double, dim> r) {
				return (-0.5 * trial_wf_2nd(r) + potential_energy(r) * trial_wf(r)) / trial_wf(r);
			};

			unsigned int n_steps(1e6);
			std::array<double, dim> new_point;
			std::vector<double> energy_samples, sequence;
			for(unsigned int n = 0; n < n_steps; ++n)
			{
				new_point = metro.step(pdf, rng);
				sequence.push_back(new_point[0]);
				energy_samples.push_back(energy(new_point));
			}

			// Calculate the final average value and its uncertainty,
			// then append those values to the output file.
			std::vector<std::vector<double>> energy_blocks(block_statistics(energy_samples, n_blocks));
            output_file << std::setw(col_width) << mu
                        << std::setw(col_width) << sigma
                        << std::setw(col_width) << energy_blocks.back()[0]
                        << std::setw(col_width) << energy_blocks.back()[1]
						<< "\n";

			calc++;
			std::cerr << "\r" << calc << " of " << std::pow(n_grid_points, 2) << " done.";
		}
	}

    output_file.close();

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
