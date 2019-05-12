#include <iostream>
#include <algorithm>
#include <fstream>
#include <cmath>
#include <string>
#include <cassert>
#include <iomanip>
#include <vector>
#include "random.hh"
#include "metropolis_NVT.hh"

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
                input_parameters_file(particle_type + ".dat");

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
    std::ofstream output(particle_type + "_output.dat");

    unsigned int n_conf(1),
				 snapshot_steps(10);

    std::vector<double> potential_en_density,
                        pressure;

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
        potential_en_density.push_back(dynamo.get_potential_energy_density());
        pressure.push_back(dynamo.get_pressure());
    }
    std::cerr << std::endl;

    dynamo.write_config("config.final");

    for(unsigned int j = 0; j < potential_en_density.size(); ++j)
    {
        output << potential_en_density[j] << " "
               << pressure[j]             << "\n";
    }
    output.close();

    std::cerr << "Acceptance rate: " << dynamo.get_acceptance_rate() << std::endl;

    return 0;
}
