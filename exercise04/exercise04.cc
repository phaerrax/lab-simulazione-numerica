#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <vector>
#include "random.hh"
#include "molecular_dynamics_sim.hh"

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

    // Initialisation procedure
    // ========================
    std::cout << "Classic Lennard-Jones fluid\n"
              << "Molecular dynamics simulation in NVE ensemble\n\n"
              << "Interatomic potential V(r) = 4 * [(1/r)^12 - (1/r)^6]\n\n"
              << "The program uses Lennard-Jones units." << std::endl;

    // Read the input parameters from a file.
    // They have to be given in a very specific way...
    std::string input_parameters_file("input.dat");
    std::cout << "Read input parameters from " << input_parameters_file << "." << std::endl;
    std::ifstream input_parameters(input_parameters_file);

    double temperature,
           n_particles,
           particle_density,
           distance_cutoff,
           time_step;
    unsigned int n_steps,
                 print_steps;

    input_parameters >> temperature
                     >> n_particles
                     >> particle_density
                     >> distance_cutoff
                     >> time_step
                     >> n_steps
                     >> print_steps;

    input_parameters.close();

    // We assume that the system is confined in a cubic cell, which volume
    // is determined by the number of particles and their density.
    // The periodic boundary conditions are already implemented in the
    // simulator class.
    double total_volume = static_cast<double>(n_particles) / particle_density;
    double cell_edge_length = std::pow(total_volume, 1. / 3.);

    molecular_dynamics_sim dynamo("config.0", cell_edge_length);
    // The coordinates in the input files are given in units of the cell
    // edge length (which is itself expressed in "Lennard-Jones units"): the
    // extra parameter in the constructor makes sure that in the simulation
    // the coordinates are properly rescaled.

    dynamo.set_particle_number(n_particles);
    dynamo.set_particle_density(particle_density);
    dynamo.set_distance_cutoff(distance_cutoff);
    dynamo.set_integration_step(time_step);

    std::cout << "Number of particles: "           << n_particles << "\n"
              << "Density of particles: "          << particle_density << "\n"
              << "Volume of the simulation cell: " << total_volume << "\n"
              << "Edge of the simulation cell: "   << cell_edge_length << "\n"
              << "The program integrates Newton equations with the Verlet method, using a time step of " << time_step << " for a total of " << n_steps << " steps." << std::endl;

    // Velocity initialisation
    // =======================
    // Generate uniformly distributed velocities, then rescale them in order
    // to match the input temperature.
    dynamo.initialise_uniform(rng);
    dynamo.rescale_velocity(temperature);

    // Integration of the equations of motion
    // ======================================
    // Now everything is in place for the algorithm to start the integration.
    // An intermediate snapshot of the system (a list of the position of the
    // particles and some thermodynamical quantities) is printed every
    // 'print_steps' only, to save some time and memory.
    std::ofstream potential_en_output("output_potential_en.dat"),
                  kinetic_en_output("output_kinetic_en.dat"),
                  temperature_output("output_temperature.dat"),
                  total_en_output("output_total_en.dat");

    double current_potential_en_density,
           current_kinetic_en_density,
           current_temperature,
           current_total_en_density;

    unsigned int n_conf = 1;
    for(unsigned int step = 1; step <= n_steps; ++step)
    {
        dynamo.move();
        if(step % print_steps == 0)
        {
            std::cout << "Number of time-steps: " << step << std::endl;

            current_potential_en_density = dynamo.get_potential_energy_density();
            current_kinetic_en_density   = dynamo.get_kinetic_energy_density();
            current_temperature          = dynamo.get_temperature();
            current_total_en_density     = current_potential_en_density + current_kinetic_en_density;

            potential_en_output << current_potential_en_density << std::endl;
            kinetic_en_output   << current_kinetic_en_density   << std::endl;
            temperature_output  << current_temperature          << std::endl;
            total_en_output     << current_total_en_density     << std::endl;

            // The integer n_conf distinguishes between the "snapshots"
            // of the system taken at regular times during the simulation.
            // The list of config_N.xyz files will be read by an external
            // program, i.e. ovito, which will be able then to visualise
            // the evolution of the system.
            dynamo.write_config_xyz("frames/config_" + std::to_string(n_conf) + ".xyz");
            ++n_conf;
        }
    }
    dynamo.write_config("config.final");

    return 0;
}
