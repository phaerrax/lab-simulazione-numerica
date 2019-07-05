#include <cmath>
#include <iomanip>
#include <map>
#include <iostream>
#include <fstream>
#include <iterator>
#include <numeric>
#include <functional>
#include <array>
#include <vector>
#include "metropolis.hh"
#include "metropolis_uniform.hh"
#include "metropolis_normal.hh"
#include "random.hh"
#include "statistics.hh"

int main(int argc, char ** argv)
{
	if(argc != 2)
	{
		std::cerr << "Error: invalid input.\n"
				  << "Syntax: " << argv[0] << " <input_file_folder>" << std::endl;
		return 1;
	}
    const unsigned int dim = 3;
	unsigned int equilibration_steps_1s,
				 equilibration_steps_2p;
    std::array<double, dim> start;

	std::string directory(argv[1]);
	directory += "/";
	std::string input_filename(directory + "input.dat");
	std::ifstream input_file(input_filename);
	if(!input_file.is_open())
	{
		std::cerr << "Unable to open " << input_filename << "." << std::endl;
		return 2;
	}
	input_file >> start[0] >> start[1] >> start[2]
			   >> equilibration_steps_1s >> equilibration_steps_2p;

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

	metropolis_uniform<dim> metro1s_uniform(-1.22, 1.22);
	metropolis_normal<dim>  metro1s_normal(0.75);
	metropolis_uniform<dim> metro2p_uniform(-2.85, 2.85);
	metropolis_normal<dim>  metro2p_normal(1.8);

	metro1s_uniform.set_starting_point(start);
	metro1s_normal.set_starting_point(start);
	metro2p_uniform.set_starting_point(start);
	metro2p_normal.set_starting_point(start);

    auto radius = [](const std::array<double, dim> & x)
    {
        return std::sqrt(
                std::inner_product(x.begin(), x.end(), x.begin(), 0.)
                );
    };

    std::function<double (std::array<double, dim>)> f1s = [radius](std::array<double, dim> x)
    {
        return std::exp(-2. * radius(x)) / M_PI;
    };
    std::function<double (std::array<double, dim>)> f2p = [radius](std::array<double, dim> x)
    {
        return std::pow(x[2], 2) * std::exp(-radius(x)) / (32 * M_PI);
    };
    // The functions to be sampled.

	// Equilibration phase
    for(unsigned int n = 0; n < equilibration_steps_1s; ++n)
    {
        metro1s_uniform.step(f1s, rng);
		metro1s_normal.step(f1s, rng);
	}

    for(unsigned int n = 0; n < equilibration_steps_2p; ++n)
	{
        metro2p_uniform.step(f2p, rng);
		metro2p_normal.step(f2p, rng);
    }

	// Actual simulation
    unsigned int n_steps(1e6);
    std::vector<std::array<double, dim>> sequence1s_uniform,
										 sequence1s_normal,
										 sequence2p_uniform,
										 sequence2p_normal;
    std::array<double, dim> new_point;
    std::vector<double> r1s_uniform,
						r1s_normal,
						r2p_uniform,
						r2p_normal;
	/*
	   Unfortunately the disk quota is too small to save every
	   sampled position, therefore I will record it only after
	   a certain number of steps.
	 */
	const unsigned int record_steps(20);
    for(unsigned int n = 0; n < n_steps; ++n)
    {
        new_point = metro1s_uniform.step(f1s, rng);
        r1s_uniform.push_back(radius(new_point));
		if(n % record_steps == 0)
			sequence1s_uniform.push_back(new_point);

		new_point = metro1s_normal.step(f1s, rng);
		r1s_normal.push_back(radius(new_point));
		if(n % record_steps == 0)
			sequence1s_normal.push_back(new_point);

        new_point = metro2p_uniform.step(f2p, rng);
        r2p_uniform.push_back(radius(new_point));
		if(n % record_steps == 0)
			sequence2p_uniform.push_back(new_point);

		new_point = metro2p_normal.step(f2p, rng);
		r2p_normal.push_back(radius(new_point));
		if(n % record_steps == 0)
			sequence2p_normal.push_back(new_point);
    }

    std::cerr << "[1s-uniform] Acceptance rate after " << n_steps << " steps: " << metro1s_uniform.get_acceptance_rate() << "." << std::endl;
    std::cerr << "[1s-normal]  Acceptance rate after " << n_steps << " steps: " << metro1s_normal.get_acceptance_rate()  << "." << std::endl;
    std::cerr << "[2p-uniform] Acceptance rate after " << n_steps << " steps: " << metro2p_uniform.get_acceptance_rate() << "." << std::endl;
    std::cerr << "[2p-normal]  Acceptance rate after " << n_steps << " steps: " << metro2p_normal.get_acceptance_rate()  << "." << std::endl;

    // Output the sequence of sampled points to a file.
    std::ofstream output_file(directory + "1s_sampled_points_uniform.dat");

    const unsigned int col_width = 16;
    output_file.precision(4);
    output_file << std::scientific;

    for(const auto & row : sequence1s_uniform)
    {
        for(auto x = row.begin(); x != row.end(); ++x)
            output_file << std::setw(col_width) << *x;
        output_file << "\n";
    }
    output_file << std::endl;
    output_file.close();

    output_file.open(directory + "2p_sampled_points_uniform.dat");
    for(const auto & row : sequence2p_uniform)
    {
        for(auto x = row.begin(); x != row.end(); ++x)
            output_file << std::setw(col_width) << *x;
        output_file << "\n";
    }
    output_file << std::endl;
    output_file.close();

    // Output the progressive values of the average radius and its standard
    // deviation, obtained with a blocking technique, to a file.
	const unsigned int block_size(30);
    std::vector<double> r1s_uniform_avg,
                        r1s_uniform_std,
                        r1s_normal_avg,
                        r1s_normal_std,
                        r2p_uniform_avg,
                        r2p_uniform_std,
                        r2p_normal_avg,
                        r2p_normal_std;

    block_statistics(
            std::begin(r1s_uniform),
            std::end(r1s_uniform),
            std::back_inserter(r1s_uniform_avg),
            std::back_inserter(r1s_uniform_std),
			block_size
            );
    block_statistics(
            std::begin(r1s_normal),
            std::end(r1s_normal),
            std::back_inserter(r1s_normal_avg),
            std::back_inserter(r1s_normal_std),
			block_size
            );
    block_statistics(
            std::begin(r2p_uniform),
            std::end(r2p_uniform),
            std::back_inserter(r2p_uniform_avg),
            std::back_inserter(r2p_uniform_std),
			block_size
            );
    block_statistics(
            std::begin(r2p_normal),
            std::end(r2p_normal),
            std::back_inserter(r2p_normal_avg),
            std::back_inserter(r2p_normal_std),
			block_size
            );

    output_file.open(directory + "1s_radius_avg.dat");
	output_file << "uniform_avg uniform_std normal_avg normal_std\n";
    for(unsigned int i = 0; i < r1s_uniform_avg.size(); ++i)
        output_file << std::setw(col_width) << r1s_uniform_avg[i]
                    << std::setw(col_width) << r1s_uniform_std[i]
                    << std::setw(col_width) << r1s_normal_avg[i]
                    << std::setw(col_width) << r1s_normal_std[i]
                    << "\n";
    output_file << std::endl;
    output_file.close();

    output_file.open(directory + "2p_radius_avg.dat");
	output_file << "uniform_avg uniform_std normal_avg normal_std\n";
    for(unsigned int i = 0; i < r2p_uniform_avg.size(); ++i)
        output_file << std::setw(col_width) << r2p_uniform_avg[i]
                    << std::setw(col_width) << r2p_uniform_std[i]
                    << std::setw(col_width) << r2p_normal_avg[i]
                    << std::setw(col_width) << r2p_normal_std[i]
                    << "\n";
    output_file << std::endl;
    output_file.close();

	rng.SaveSeed();
    return 0;
}
