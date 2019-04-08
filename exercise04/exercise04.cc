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
           time_step,
           n_steps,
           print_steps;

    input_parameters >> temperature
                     >> n_particles
                     >> particle_density
                     >> distance_cutoff
                     >> time_step
                     >> n_steps
                     >> print_steps;

    input_parameters.close();

    molecular_dynamics_sim dynamo("config.0");

    dynamo.set_particle_number(n_particles);
    dynamo.set_particle_density(particle_density);
    dynamo.set_distance_cutoff(distance_cutoff);
    dynamo.set_integration_step(time_step);

    total_volume = static_cast<double>(n_particles) / particle_density;
    cell_edge_length = std::pow(total_volume, 1. / 3.);

    std::cout << "Number of particles: "           << n_particles << "\n"
              << "Density of particles: "          << particle_density << "\n"
              << "Volume of the simulation cell: " << total_volume << "\n"
              << "Edge of the simulation cell: "   << cell_edge_length << "\n"
              << "The program integrates Newton equations with the Verlet method, using a time step of " << time_step << " for a total of " << n_steps << " steps." << std::endl;

    // Velocity initialisation: generate uniformly distributed velocities,
    // then rescale them to match a desired temperature.
    dynamo.initialise_uniform(rng);
    dynamo.rescale_velocity(temperature);

    int nconf = 1;
    for(int istep=1; istep <= nstep; ++istep)
    {
        Move();           //Move particles with Verlet algorithm
        if(istep%iprint == 0)
            std::cout << "Number of time-steps: " << istep << std::endl;
        if(istep%10 == 0)
        {
            Measure();     //Properties measurement
            //        ConfXYZ(nconf);//Write actual configuration in XYZ format //Commented to avoid "filesystem full"! 
            nconf += 1;
        }
    }
    ConfFinal();         //Write final configuration to restart

    return 0;
}
