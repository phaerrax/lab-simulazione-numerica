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
    Random rng1, rng2;
    int seed[4];
    int p1, p2, p3, p4;
    std::ifstream Primes("Primes");
    if(Primes.is_open())
    {
        Primes >> p1 >> p2
               >> p3 >> p4;
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
                rng1.SetRandom(seed,p1,p2);
                rng2.SetRandom(seed,p3,p4);
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
              << "The program uses Lennard-Jones units.\n"
              << "-----------------------------------------------------\n"
              << "The program will launch two separate simulations,\n"
              << "using different istantiations of the random number\n"
              << "generator, in order to be able to analyse the warm-up\n"
              << "period of the Markov chain." << std::endl;

    // Read the input parameters from a file.
    // They have to be given in a very specific way...
    std::cout << "Read input parameters from " << input_parameters_file << "." << std::endl;
    std::ifstream input_parameters(input_parameters_file);

    double temperature,
           n_particles,
           particle_density,
           distance_cutoff,
		   step_stdev;
    unsigned int n_steps;

    input_parameters >> temperature
                     >> n_particles
                     >> particle_density
                     >> distance_cutoff
					 >> step_stdev
                     >> n_steps;

    input_parameters.close();

    // We assume that the system is confined in a cubic cell, which volume
    // is determined by the number of particles and their density.
    // The periodic boundary conditions are already implemented in the
    // simulator class.
    double total_volume = static_cast<double>(n_particles) / particle_density;
    double cell_edge_length = std::pow(total_volume, 1. / 3.);

    metropolis_NVT dynamo1("config.0", cell_edge_length),
                   dynamo2("config.fcc", cell_edge_length);
    // The coordinates in the input files are given in units of the cell
    // edge length (which is itself expressed in "Lennard-Jones units"): the
    // extra parameter in the constructor makes sure that in the simulation
    // the coordinates are properly rescaled.

    // Initialisation
    // ==============
    dynamo1.set_particle_density(particle_density);
    dynamo1.set_step_stdev(step_stdev);
    dynamo1.set_temperature(temperature);
    dynamo1.set_distance_cutoff(distance_cutoff);
    dynamo2.set_particle_density(particle_density);
    dynamo2.set_step_stdev(step_stdev);
    dynamo2.set_temperature(temperature);
    dynamo2.set_distance_cutoff(distance_cutoff);

    std::cout << "Number of particles: "           << n_particles << "\n"
              << "Density of particles: "          << particle_density << "\n"
              << "Temperature: "                   << temperature << "\n"
              << "Volume of the simulation cell: " << total_volume << "\n"
              << "Edge of the simulation cell: "   << cell_edge_length << std::endl;

    double progress;

    // Simulation
    // ==========
    std::ofstream output("warmup_analysis_output.dat");

    std::vector<double> potential_en_density1,
                        potential_en_density2;

    for(unsigned int step = 1; step <= n_steps; ++step)
    {
		// A "Monte Carlo step" consists in moving all the particles.
		for(unsigned int i = 0; i < n_particles; ++i)
        {
			dynamo1.next(rng1);
			dynamo2.next(rng2);
        }

        if(step % 1000)
        {
            progress = 100 * static_cast<double>(step) / n_steps;
            std::cerr << "\rNumber of time-steps: " << step << " / " << n_steps << " (" << std::round(progress) << "%)";
        }

        progress = 100 * static_cast<double>(step) / n_steps;
        potential_en_density1.push_back(dynamo1.get_potential_energy_density());
        potential_en_density2.push_back(dynamo2.get_potential_energy_density());
    }
    std::cerr << std::endl;

    for(unsigned int j = 0; j < potential_en_density1.size(); ++j)
    {
        output << potential_en_density1[j] << " "
               << potential_en_density2[j] << "\n";
    }
    output.close();

    std::cerr << "Acceptance rate: " << dynamo1.get_acceptance_rate() 
              << "and " << dynamo2.get_acceptance_rate() << std::endl;

    return 0;
}
