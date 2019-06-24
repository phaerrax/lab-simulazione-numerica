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
#include "statistics.hh"

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
    if(argc != 3)
    {
        std::cerr << "Error: invalid input.\n"
				  << "Syntax: " << argv[0] << " <element> <phase>\n"
				  << "to use the parameters in the file \"./<element>/<phase>/input.dat\"." << std::endl;
        std::cerr << "Error: too few arguments." << std::endl;
        return 1;
    }

	std::string particle_type(argv[1]),
				phase(argv[2]);

    // The directory in which the input parameters and output files are
    // stored.
    std::string prefix = particle_type + "/" + phase + "/";

    std::cout << "Classic Lennard-Jones fluid: "
              << particle_type << " in a " << phase << " phase.\n"
              << "Molecular dynamics simulation in NVT ensemble\n\n"
              << "Interatomic potential V(r) = 4 * [(1/r)^12 - (1/r)^6]\n\n"
              << "The program uses Lennard-Jones units." << std::endl;

    // Read the input parameters from a file.
    // They have to be given in a very specific way...
    std::string input_parameters_file(prefix + "input.dat");
    std::cout << "Read input parameters from " << input_parameters_file << "." << std::endl;
    std::ifstream input_parameters(input_parameters_file);

    double temperature,
           n_particles,
           particle_density,
           distance_cutoff,
		   step_stdev;
    unsigned int n_steps,
				 n_equilibration_steps,
				 block_size;

    input_parameters >> temperature
                     >> n_particles
                     >> particle_density
                     >> distance_cutoff
					 >> step_stdev
                     >> n_steps
					 >> block_size
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
    unsigned int n_conf(1),
				 snapshot_steps(10),
	             n_bins(100); // of the radial dist function.

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
		/*
		   I have to disable this part since the files it generates
		   quickly eat up all the available disk space.
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
		*/

		auto results = dynamo.measure(n_bins, max_distance_rd);
        potential_en_density.push_back(std::get<0>(results));
        pressure.push_back(std::get<1>(results));
		radial_distribution.push_back(std::move(std::get<2>(results)));
		// The results tuple is deleted anyway at the end of the iteration.
    }
    std::cerr << std::endl;

    dynamo.write_config("config.final");

    // Calculate average values and standard deviation of the measured
    // physical quantities.
    std::vector<double> potential_en_density_avg,
                        potential_en_density_std,
                        pressure_avg,
                        pressure_std;

    block_statistics(
            std::begin(potential_en_density),
            std::end(potential_en_density),
            std::back_inserter(potential_en_density_avg),
            std::back_inserter(potential_en_density_std),
            block_size
            );

    block_statistics(
            std::begin(pressure),
            std::end(pressure),
            std::back_inserter(pressure_avg),
            std::back_inserter(pressure_std),
            block_size
            );

	const std::string suffix("_metropolis.dat");
	std::ofstream avg_rd_output(prefix + "avg_radial_dist" + suffix),
	              final_rd_output(prefix + "histogram_radial_dist" + suffix),
                  pot_en_output(prefix + "pot_energy" + suffix),
                  pressure_output(prefix + "pressure" + suffix);
	pot_en_output      << "steps pot_en_avg pot_en_std\n";
	pressure_output    << "steps pressure_avg pressure_std\n";
	unsigned int steps;
    for(unsigned int j = 0; j < potential_en_density_avg.size(); ++j)
    {
		// All four vectors in this loop have the same length,
		// that is the number of blocks.
        steps = (j + 1) * block_size;
		pot_en_output      << steps                       << " "
						   << potential_en_density_avg[j] << " "
					       << potential_en_density_std[j] << "\n";
        pressure_output    << steps                       << " "
			               << pressure_avg[j]             << " "
                           << pressure_std[j]             << "\n";
    }
    pot_en_output.close();
    pressure_output.close();

	// Calculate the average radial distribution with a blocking technique.
	// Instead of summing directly each vector...
	// ... transpose the "matrix", obtaining a vector of subvectors, each
	// subvector containing the values, step by step, assumed in a bin.
	// Then, with the results of the average and the standard deviation,
	// build a vector back.
	std::vector<std::vector<double>> tr_radial_distribution(
			radial_distribution[0].size()
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
	// that is on each "column" of radial_distribution.
	std::vector<double> rd_final_stdev;
	std::vector<std::vector<double>> rd_averages;

	for(auto & row : tr_radial_distribution)
	{
        std::vector<double> avg, err;
		block_statistics(
                std::begin(row),
                std::end(row),
                std::back_inserter(avg),
                std::back_inserter(err),
				block_size
                );
		// I am only interested in the averages and in the final stdev.

        // avg will be a n_blocks-long vector; its elements
        // are the block averages of the i-th bin.
        // rd_averages itself has n_bins elements.
		rd_averages.push_back(std::move(avg));
		rd_final_stdev.push_back(err.back());
	}

	// Output the block averages of the radial distribution function.
	// The first row contains the upper bound of each bin...
	for(auto elem : bin_upper_bounds)
		avg_rd_output << elem << " ";
	avg_rd_output << "\n";
	// ...then each row that follows contains the average, block by
	// block, of each bin.

	/*
	   - rd_averages[i][*] are the values of the i-th bin
	   as the block number increases.
	   - rd_averages[*][j] are the values of the bins
	   as the j-th data block is added to the count.
	   - rd_averages[j][i] is the average value of the j-th
	   bin at the i-th data block.

	   rd_averages.size() == n_bins,
	   rd_averages[0].size() == n_blocks (whatever that is).

	   I want to output the values in such a way that each
	   line corresponds to a block, i.e. it contains the
	   values of all bins, sequentially, at a particular
	   block.
	   The following line will contain again the value of 
	   all bins, but evaluated with one data block more
	   (therefore they are a more accurate measurement).

	   For this reason, I need to cycle on rd_averages
	   firstly on the block index, secondly on the bin 
	   index.
	 */
	for(unsigned int j = 0; j < rd_averages[0].size(); ++j)
	{
		for(unsigned int i = 0; i < rd_averages.size(); ++i)
			avg_rd_output << rd_averages[i][j] << "  ";
		avg_rd_output << "\n";
	}
	avg_rd_output.close();

	// Output the "final" radial distribution function, with its uncertainty.
	// The first column contains the upper bound of each bin.
	final_rd_output << "bin_upper_bound final_rd_avg final_rd_std\n";
    for(unsigned int j = 0; j < rd_final_stdev.size(); ++j)
	{
		final_rd_output << bin_upper_bounds[j] << " "
			            << rd_averages[j].back() << " "
		  	            << rd_final_stdev[j] << "\n";
	}
	final_rd_output.close();

    std::cerr << "Acceptance rate: " << dynamo.get_acceptance_rate() << std::endl;

	rng.SaveSeed();
    return 0;
}
