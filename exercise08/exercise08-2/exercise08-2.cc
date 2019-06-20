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
#include "statistics.hh"

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

	if(argc != 3)
	{
		std::cerr << "Too few arguments (3 required)." << std::endl;
		return 1;
	}

    const unsigned int dim = 1;
	double mu(std::stod(argv[1])),
		   sigma(std::stod(argv[2]));
    double L(1.2);
    metropolis_uniform<dim> metro(-L, L);
    
    std::array<double, dim> start = {0};
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
		double plus  = std::exp(-0.5 * std::pow((x + mu) / sigma, 2)),
			   minus = std::exp(-0.5 * std::pow((x - mu) / sigma, 2));
		return std::pow(sigma, -2) * (
				-(plus + minus) +
				std::pow(sigma, -2) * (std::pow(x - mu, 2) * minus + std::pow(x + mu, 2) * plus)
				);
	};

	auto potential_energy = [](std::array<double, dim> r) {
		double x = r[0];
		return (std::pow(x, 2) - 2.5) * std::pow(x, 2);
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

    std::cerr << "Acceptance rate after " << n_steps << " steps: " << metro.get_acceptance_rate() << "." << std::endl;

    // Output the progressive values of the average energy and its standard
    // deviation, obtained with a blocking technique, to a file.
	std::ofstream output_file("energy.dat");
    const unsigned int col_width = 16;
    output_file.precision(4);
    output_file << std::scientific;

	unsigned int steps,
				 n_blocks(100);
    std::vector<double> energy_avg,
                        energy_std;
    block_statistics(
            std::begin(energy_samples),
            std::end(energy_samples),
            std::back_inserter(energy_avg),
            std::back_inserter(energy_std),
            energy_samples.size() / n_blocks
            );
    for(unsigned int j = 0; j < energy_samples.size(); ++j)
	{
		steps = (j + 1) * n_steps / n_blocks;
        output_file << std::setw(col_width) << steps
                    << std::setw(col_width) << energy_avg[j]
                    << std::setw(col_width) << energy_std[j]
                    << "\n";
	}
    output_file.close();

	std::ofstream sequence_output_file("drawn_points.dat");
	std::ostream_iterator<double> sequence_output(sequence_output_file, "\n");
	std::copy(sequence.begin(), sequence.end(), sequence_output);

    return 0;
}
