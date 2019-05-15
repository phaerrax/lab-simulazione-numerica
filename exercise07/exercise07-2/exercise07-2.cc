#include <iostream>
#include <algorithm>
#include <fstream>
#include <cmath>
#include <utility>
#include <string>
#include <tuple>
#include <cassert>
#include <iomanip>
#include <vector>
#include "random.hh"
#include "metropolis_NVT.hh"

std::vector<std::vector<double>> block_statistics(const std::vector<double> &, unsigned int);

int main(int argc, char *argv[])
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
    // Get the name of the input file from the command line
    if(argc != 2)
    {
        std::cerr << "Error: too few arguments." << std::endl;
        return 1;
    }
	std::string particle_type(argv[1]),
                input_parameters_file("input." + particle_type);

    std::cout << "Classic Lennard-Jones fluid: "
              << particle_type << "\n"
              << "Molecular dynamics simulation in NVT ensemble\n\n"
              << "Interatomic potential V(r) = 4 * [(1/r)^12 - (1/r)^6]\n\n"
              << "The program uses Lennard-Jones units." << std::endl;

    // Read the input parameters from a file.
    // They have to be given in a very specific way...
    std::cout << "Read input parameters from " << input_parameters_file << "." << std::endl;
    std::ifstream input_parameters(input_parameters_file);

    double temperature,
           n_particles,
           particle_density,
           distance_cutoff,
		   step_stdev;
    unsigned int n_steps,
				 n_equilibration_steps;

    input_parameters >> temperature
                     >> n_particles
                     >> particle_density
                     >> distance_cutoff
					 >> step_stdev
                     >> n_steps
					 >> n_equilibration_steps;

    input_parameters.close();

    // We assume that the system is confined in a cubic cell, which volume
    // is determined by the number of particles and their density.
    // The periodic boundary conditions are already implemented in the
    // simulator class.
    double total_volume = static_cast<double>(n_particles) / particle_density;
    double cell_edge_length = std::pow(total_volume, 1. / 3.);

    metropolis_NVT dynamo("config.0", cell_edge_length);
    // The coordinates in the input files are given in units of the cell
    // edge length (which is itself expressed in "Lennard-Jones units"): the
    // extra parameter in the constructor makes sure that in the simulation
    // the coordinates are properly rescaled.

    // Initialisation
    // ==============
    dynamo.set_particle_density(particle_density);
    dynamo.set_step_stdev(step_stdev);
    dynamo.set_temperature(temperature);
    dynamo.set_distance_cutoff(distance_cutoff);

    std::cout << "Number of particles: "           << n_particles << "\n"
              << "Density of particles: "          << particle_density << "\n"
              << "Temperature: "                   << temperature << "\n"
              << "Volume of the simulation cell: " << total_volume << "\n"
              << "Edge of the simulation cell: "   << cell_edge_length << std::endl;

    // Equilibration run
    // =================
    double progress;
    std::cerr << "Equilibrating...";
    for(unsigned int step = 0; step < n_equilibration_steps; ++step)
    {
		// A "Monte Carlo step" consists in moving all the particles.
		for(unsigned int i = 0; i < n_particles; ++i)
			dynamo.next(rng);

        if(step % 100 == 0)
        {
            progress = 100 * static_cast<double>(step) / n_steps;
            std::cerr << "\rEquilibrating... " << std::round(progress) << "%";
        }
    }
    std::cerr << "\rEquilibrating... done." << std::endl;

    // Simulation
    // ==========
    //std::ofstream output(particle_type + "_output.dat");
	std::ofstream avg_rd_output(particle_type + "_output.gofr.0"),
	              final_rd_output(particle_type + "_output.gave.0"),
	              other_output(particle_type + "_output.0");

    unsigned int n_conf(1),
				 snapshot_steps(10);

	unsigned int n_bins(100); // of the radial dist function.

    //std::vector<double> potential_en_density,
    //                    pressure;

	std::vector<std::vector<double>> radial_distribution;
	std::vector<double> potential_en_density,
						pressure;

	radial_distribution.reserve(n_steps);
	potential_en_density.reserve(n_steps);
	pressure.reserve(n_steps);

	double max_distance_rd(0.5 * cell_edge_length);
	// Stop at this distance when calculating the radial distribution.
	std::vector<double> bin_upper_bounds(n_bins);
	for(unsigned int i = 0; i < n_bins; ++i)
		bin_upper_bounds[i] = max_distance_rd / n_bins * (i + 1);

	// Write the first (equilibrated) configuration on a file.
    dynamo.write_config_xyz("frames/config_0.xyz");

    for(unsigned int step = 1; step <= n_steps; ++step)
    {
		// A "Monte Carlo step" consists in moving all the particles.
		for(unsigned int i = 0; i < n_particles; ++i)
			dynamo.next(rng);

        if(step % 1000)
        {
            progress = 100 * static_cast<double>(step) / n_steps;
            std::cerr << "\rNumber of time-steps: " << step << " / " << n_steps << " (" << std::round(progress) << "%)";
        }

        progress = 100 * static_cast<double>(step) / n_steps;
        if(step % snapshot_steps == 0)
        {
            // The integer n_conf distinguishes between the "snapshots"
            // of the system taken at regular times during the simulation.
            // The list of config_N.xyz files will be read by an external
            // program, i.e. ovito, which will be able then to visualise
            // the evolution of the system.
            dynamo.write_config_xyz("frames/config_" + std::to_string(n_conf) + ".xyz");
            ++n_conf;
        }

		auto results = dynamo.measure(n_bins, max_distance_rd);
        potential_en_density.push_back(std::get<0>(results));
        pressure.push_back(std::get<1>(results));
		radial_distribution.push_back(std::move(std::get<2>(results)));
		// The results tuple is deleted anyway at the end of the iteration.
    }
    std::cerr << std::endl;

    dynamo.write_config("config.final");

	// Calculate the average radial distribution with a blocking technique.
	// Instead of summing directly each vector...
	// ... transpose the "matrix", obtaining a vector of subvectors, each
	// subvector containing the values, step by step, assumed in a bin.
	// Then, with the results of the average and the standard deviation,
	// build a vector back.
	std::vector<std::vector<double>> tr_radial_distribution(
			radial_distribution.begin()->size()
			);
	// radial_distribution.size() == n_steps,
	// radial_distribution[0].size() == n_bins.
	for(auto & row : tr_radial_distribution)
		row.resize(radial_distribution.size());
	// tr_radial_distribution.size() == n_bins,
	// tr_radial_distribution[0].size() == n_steps.

	for(unsigned int i = 0; i < tr_radial_distribution.size(); ++i)
		for(unsigned int j = 0; j < radial_distribution.size(); ++j)
			tr_radial_distribution[i][j] = radial_distribution[j][i];

	// Call block_statistics on each subvector of tr_radial_distribution,
	// that is on each "column" or radial_distribution.
	std::vector<double> rd_final_stdev(tr_radial_distribution.size());
	std::vector<std::vector<double>> rd_averages(tr_radial_distribution.size());

	for(unsigned int i = 0; i < tr_radial_distribution.size(); ++i)
	{
		auto v = block_statistics(tr_radial_distribution[i], 100);
		// For each bin of the histogram, we get a vector of pairs
		// (average, stdev). I am only interested in the averages
		// and in the final stdev.
		for(const auto & pair : v)
			rd_averages[i].push_back(pair.front());
			// rd_averages[i] will be a n_blocks-long vector; its elements
			// are the block averages of the i-th bin.
			// rd_averages itself has n_bins elements.
		rd_final_stdev[i] = v.back().back();
	}

	// Output the block averages of the radial distribution function.
	// The first row contains the upper bound of each bin...
	for(auto elem : bin_upper_bounds)
		avg_rd_output << elem << " ";
	avg_rd_output << "\n";
	// ...then each row that follows contains the average, block by
	// block, of each bin.
	for(unsigned int i = 0; i < rd_averages.size(); ++i)
	{
		// Each element in rd_averages is a histogram, that is the average
		// histogram for that block.
		// rd_averages[j][i] is the average value of the j-th bin for the
		// i-th block.
		for(unsigned int j = 0; j < rd_averages[0].size(); ++j)
			// Each element in row is the value in a bin.
			// i = 0, ..., n_bins;
			// j = 0, ..., n_blocks.
			avg_rd_output << rd_averages[j][i] << "  ";
		avg_rd_output << "\n";
	}
	avg_rd_output.close();

	// Output the "final" radial distribution function, with its uncertainty.
	// The first column contains the upper bound of each bin.
	final_rd_output << "bin_upper_bound final_rd_avg final_rd_std\n";
    for(unsigned int j = 0; j < rd_final_stdev.size(); ++j)
		final_rd_output << bin_upper_bounds[j] << " "
			            << rd_averages[j].back() << " "
		  	            << rd_final_stdev[j] << "\n";
	final_rd_output.close();

	other_output << "potential_en_density pressure\n";
    for(unsigned int j = 0; j < potential_en_density.size(); ++j)
    {
        other_output << potential_en_density[j] << " "
                     << pressure[j]             << "\n";
    }
    other_output.close();

    std::cerr << "Acceptance rate: " << dynamo.get_acceptance_rate() << std::endl;

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
