#include <iostream>
#include <algorithm>
#include <fstream>
#include <cmath>
#include <string>
#include <cassert>
#include <vector>
#include "random.hh"
#include "statistics.hh"
#include "molecular_dynamics_sim.hh"

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
    // Gather the type of the particle in the system and its thermodinamical
    // phase from the command-line arguments.
    if(argc != 2)
    {
        std::cerr << "Error: invalid input.\n\n"
				  << "Syntax: " << argv[0] << " <phase>\n"
				  << "(for example: " << argv[0] << " liquid)"
				  << std::endl;
        return 1;
    }

    std::string phase(argv[1]);

    // The directory in which the input parameters and output files are
    // stored.
    std::string prefix = phase + "/";

    std::cout << "Classic Lennard-Jones fluid in a " << phase << " phase.\n"
              << "Molecular dynamics simulation in NVE ensemble\n\n"
              << "Interatomic potential V(r) = 4 * [(1/r)^12 - (1/r)^6]\n\n"
              << "The program uses Lennard-Jones units." << std::endl;

    // Read the input parameters from a file.
    // They have to be given in a very specific way...
    std::string input_parameters_file(prefix + "input.dat");
    std::cout << "Reading input parameters from " << input_parameters_file << "..." << std::endl;
    std::ifstream input_parameters(input_parameters_file);

    double input_temperature,
           n_particles,
           particle_density,
           distance_cutoff,
           time_step;
    unsigned int n_steps,
                 snapshots_steps,
                 measurements_steps,
				 equilibration_steps,
                 block_size;

	if(input_parameters.is_open())
	{
		input_parameters >> input_temperature
						 >> n_particles
						 >> particle_density
						 >> distance_cutoff
						 >> time_step
						 >> n_steps
						 >> snapshots_steps
						 >> measurements_steps
						 >> equilibration_steps
						 >> block_size;
	}
	else
	{
		std::cerr << "Error: unable to open " << input_parameters_file << "." << std::endl;
		return 2;
	}
    input_parameters.close();
    std::cout << "Input complete." << std::endl;

    // We assume that the system is confined in a cubic cell, which volume
    // is determined by the number of particles and their density.
    // The periodic boundary conditions are already implemented in the
    // simulator class.
    double total_volume = static_cast<double>(n_particles) / particle_density;
    double cell_edge_length = std::pow(total_volume, 1. / 3.);

    molecular_dynamics_sim dynamo_equilibration("config.0", cell_edge_length);
    // The coordinates in the input files are given in units of the cell
    // edge length (which is itself expressed in "Lennard-Jones units"): the
    // extra parameter in the constructor makes sure that in the simulation
    // the coordinates are properly rescaled.

    dynamo_equilibration.set_particle_number(n_particles);
    dynamo_equilibration.set_particle_density(particle_density);
    dynamo_equilibration.set_distance_cutoff(distance_cutoff);
    dynamo_equilibration.set_integration_step(time_step);

    std::cout << "Number of particles: "           << n_particles << "\n"
              << "Density of particles: "          << particle_density << "\n"
              << "Volume of the simulation cell: " << total_volume << "\n"
              << "Edge of the simulation cell: "   << cell_edge_length << "\n"
              << "The program integrates Newton equations with the Verlet method, using a time step of "
			  << time_step << " for a total of " << n_steps << " steps, "
			  << "after an equilibraton period of " << equilibration_steps << " steps." << std::endl;

    // Velocity initialisation
    // =======================
    // Generate uniformly distributed velocities, then rescale them in order
    // to match the input temperature.
    dynamo_equilibration.initialise_maxwellboltzmann(input_temperature, rng);
    //dynamo_equilibration.initialise_uniform(rng);
    dynamo_equilibration.rescale_velocity(input_temperature);

    // Equilibration run
    // =================
    double progress;
    std::cerr << "Equilibrating...";
    for(unsigned int step = 1; step < equilibration_steps; ++step)
    {
        dynamo_equilibration.move();
        if(step % 100 == 0)
        {
            progress = 100 * static_cast<double>(step) / equilibration_steps;
            std::cerr << "\rEquilibrating... " << std::round(progress) << "%";
        }
    }
    std::cerr << "\rEquilibrating... done." << std::endl;
    // Stop a step before the end: we need to output the second-to-last
    // configuration, which will be used later to restart the simulation.
    dynamo_equilibration.write_config("config.prefinal");

    // Move again, for the last time.
    dynamo_equilibration.move();
    dynamo_equilibration.write_config("config.final");


    // Initialise a new simulator with the final configuration of the
    // equilibration run.
    molecular_dynamics_sim dynamo("config.final", "config.prefinal", time_step, cell_edge_length);
	/*
	   From the two config files, the constructor calculates, using the
	   Verlet algorithm, the current velocity. From that, the temperature
	   can be calculated and then the velocities can be rescaled so that
	   the temperature matches the input one.
	*/

    dynamo.set_particle_number(n_particles);
    dynamo.set_particle_density(particle_density);
    dynamo.set_distance_cutoff(distance_cutoff);

    dynamo.rescale_velocity(input_temperature);

    // Integration of the equations of motion
    // ======================================
    // Now everything is in place for the algorithm to start the integration.
    // An intermediate snapshot of the system (a list of the position of the
    // particles and some thermodynamical quantities) is printed every
    // 'snapshots_steps' only, to save some time and memory.
    std::ofstream potential_en_output(prefix + "pot_energy.dat"),
                  kinetic_en_output(prefix + "kin_energy.dat"),
                  total_en_output(prefix + "tot_energy.dat"),
                  temperature_output(prefix + "temperature.dat"),
                  pressure_output(prefix + "pressure.dat");

    unsigned int n_conf(1);

    dynamo.write_config_xyz(prefix + "frames/config_0.xyz");

    std::vector<double> potential_en_density,
                        kinetic_en_density,
                        temperature,
                        pressure,
                        total_en_density;

    for(unsigned int step = 1; step <= n_steps; ++step)
    {
        dynamo.move();
        if(step % measurements_steps == 0)
        {
            progress = 100 * static_cast<double>(step) / n_steps;
            std::cerr << "\rNumber of time-steps: " << step << " / " << n_steps << " (" << std::round(progress) << "%)";

			auto results = dynamo.measure();

            temperature.push_back(std::get<0>(results));
            potential_en_density.push_back(std::get<1>(results));
            kinetic_en_density.push_back(std::get<2>(results));
            pressure.push_back(std::get<3>(results));
			// Sum potential and kinetic energy together.
			total_en_density.push_back(std::get<1>(results) + std::get<2>(results));
        }
		if(step % snapshots_steps == 0)
		{
            dynamo.write_config_xyz(prefix + "frames/config_" + std::to_string(n_conf) + ".xyz");
            ++n_conf;
		}
    }
    std::cerr << std::endl;

    dynamo.write_config("config.final");

    // Calculate average values and standard deviation of the measured
    // physical quantities.
    std::vector<double> potential_en_density_avg,
                        potential_en_density_std,
                        kinetic_en_density_avg,
                        kinetic_en_density_std,
                        total_en_density_avg,
                        total_en_density_std,
                        temperature_avg,
                        temperature_std,
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
            std::begin(kinetic_en_density),
            std::end(kinetic_en_density),
            std::back_inserter(kinetic_en_density_avg),
            std::back_inserter(kinetic_en_density_std),
            block_size
            );

    block_statistics(
            std::begin(total_en_density),
            std::end(total_en_density),
            std::back_inserter(total_en_density_avg),
            std::back_inserter(total_en_density_std),
            block_size
            );

    block_statistics(
            std::begin(temperature),
            std::end(temperature),
            std::back_inserter(temperature_avg),
            std::back_inserter(temperature_std),
            block_size
            );

    block_statistics(
            std::begin(pressure),
            std::end(pressure),
            std::back_inserter(pressure_avg),
            std::back_inserter(pressure_std),
            block_size
            );

	potential_en_output << "steps pot_en_avg pot_en_std\n";
	kinetic_en_output   << "steps kin_en_avg kin_en_std\n";
	total_en_output     << "steps tot_en_avg tot_en_std\n";
	temperature_output  << "steps temperature_avg temperature_std\n";
	pressure_output     << "steps pressure_avg pressure_std\n";
	unsigned int steps;
    for(unsigned int j = 0; j < potential_en_density_avg.size(); ++j)
    {
        steps = (j + 1) * block_size;
        potential_en_output << steps                       << " "
			                << potential_en_density_avg[j] << " "
                            << potential_en_density_std[j] << "\n";
        kinetic_en_output   << steps                       << " "
			                << kinetic_en_density_avg[j]   << " "
                            << kinetic_en_density_std[j]   << "\n";
        total_en_output     << steps                       << " "
			                << total_en_density_avg[j]     << " "
                            << total_en_density_std[j]     << "\n";
        temperature_output  << steps                       << " "
			                << temperature_avg[j]          << " "
                            << temperature_std[j]          << "\n";
        pressure_output     << steps                       << " "
			                << pressure_avg[j]             << " "
                            << pressure_std[j]             << "\n";
    }
    potential_en_output.close();
    kinetic_en_output.close();
    total_en_output.close();
    temperature_output.close();
    pressure_output.close();

    return 0;
}
