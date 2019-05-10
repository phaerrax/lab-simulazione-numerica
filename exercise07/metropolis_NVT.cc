#include "metropolis_NVT.hh"
#include <algorithm>
#include <cmath>
#include <numeric>
#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>
#include <sstream>

metropolis_NVT::metropolis_NVT(const std::string & initial_config_file, double input_cell_edge_length)
{   
    cell_edge_length = input_cell_edge_length;
    n_particles = 0;
    accepted_proposals = 0;
    total_proposals = 0;

    // Clear the previous configuration (if there is any).
    current_configuration.clear();

    // Read the initial configuration from a file.
    std::cout << "Read initial configuration from file " << initial_config_file << "." << std::endl;
    std::ifstream input_config(initial_config_file);

    std::string line, coordinate;

    // Read the file line by line.
    while(std::getline(input_config, line))
    {
        // Save each line in a new vector.
        // (The coordinates are given in units of the cell edge length, so
        // we need to rescale them.)
        std::istringstream iss(line);
		// We expect that a line consists in all coordinates of a particle.
		// The number of lines in the initial config file will be the
		// number of particles.
        while(iss >> coordinate)
            current_configuration.push_back(cell_edge_length * std::stod(coordinate));

        ++n_particles;
    }

    input_config.close();
    return;
}

void metropolis_NVT::set_temperature(double input_temperature)
{
    temperature = input_temperature;
    return;
}

void metropolis_NVT::set_particle_density(double input_particle_density)
{
    particle_density = input_particle_density;
    return;
}

void metropolis_NVT::set_distance_cutoff(double input_distance_cutoff)
{
    distance_cutoff = input_distance_cutoff;
    return;
}

void metropolis_NVT::set_step_stdev(double input_stdev)
{
    stdev = input_stdev;
    return;
}

double metropolis_NVT::size() const
{
    return current_configuration.size();
}

void metropolis_NVT::next(Random & rng)
{
	// Select a particle at random.
	unsigned int selected_particle = static_cast<unsigned int>(rng.Rannyu(0, n_particles));
    double acceptance_threshold;

	std::vector<double> selected_position {
		current_configuration[3 * selected_particle],
		current_configuration[3 * selected_particle + 1],
		current_configuration[3 * selected_particle + 2]
	};

	double x, y, z, d;
	double old_energy(0);
	for(unsigned int n = 0; n < n_particles; ++n)
		if(n != selected_particle)
		{
			x = quotient(selected_position[0] - current_configuration[3 * n]);
			y = quotient(selected_position[1] - current_configuration[3 * n + 1]);
			z = quotient(selected_position[2] - current_configuration[3 * n + 2]);
			d = std::hypot(std::hypot(x, y), z);

			if(d < distance_cutoff)
				old_energy += std::pow(d, -12) - std::pow(d, -6);
		}
	old_energy *= 4;

	// Generate a proposal. 
	std::vector<double> proposed_position {
		quotient(selected_position[0] + rng.Gauss(0, stdev)),
		quotient(selected_position[1] + rng.Gauss(0, stdev)),
		quotient(selected_position[2] + rng.Gauss(0, stdev))
	};

	// Calculate the interaction energy in the new configuration.
	double new_energy(0);
	for(unsigned int n = 0; n < n_particles; ++n)
		if(n != selected_particle)
		{
			x = quotient(proposed_position[0] - current_configuration[3 * n]);
			y = quotient(proposed_position[1] - current_configuration[3 * n + 1]);
			z = quotient(proposed_position[2] - current_configuration[3 * n + 2]);
			d = std::hypot(std::hypot(x, y), z);

			if(d < distance_cutoff)
				new_energy += std::pow(d, -12) - std::pow(d, -6);
		}

	new_energy *= 4;

	double energy_diff = new_energy - old_energy;
    total_proposals++;
	if(energy_diff < 0)
	{
		accepted_proposals++;
		current_configuration[3 * selected_particle]     = proposed_position[0];
		current_configuration[3 * selected_particle + 1] = proposed_position[1];
		current_configuration[3 * selected_particle + 2] = proposed_position[2];
	}
	else
	{
		acceptance_threshold = std::exp(-energy_diff / temperature);
		if(rng.Rannyu() < acceptance_threshold)
		{
			accepted_proposals++;
			current_configuration[3 * selected_particle]     = proposed_position[0];
			current_configuration[3 * selected_particle + 1] = proposed_position[1];
			current_configuration[3 * selected_particle + 2] = proposed_position[2];
		}
	}
}

double metropolis_NVT::get_potential_energy_density() const
{
    unsigned int n_coordinates = 3;
    std::vector<double> displacement(n_coordinates);
    double potential_en(0), sq_distance;
    // Cycle over pairs of particles only.
    for(unsigned int i = 0; i < n_particles - 1; ++i)
        for(unsigned int j = i + 1; j < n_particles; ++j)
        {
            for(unsigned int d = 0; d < n_coordinates; ++d)
                displacement[d] = quotient(current_configuration[n_coordinates * i + d] - current_configuration[n_coordinates * j + d]);
            // In the computation of the potential en include only
            // the pairs whose distance is less than the cutoff radius.
            sq_distance = std::inner_product(displacement.begin(), displacement.end(), displacement.begin(), 0.);
            // In the expression for the potential energy the distance appears
            // in even powers only, so we don't need to compute the square
            // root.
            if(sq_distance < std::pow(distance_cutoff, 2))
                potential_en += 4. * pow(sq_distance, -6) - 4. * pow(sq_distance, -3);
        }          

    return potential_en / n_particles;
}

double metropolis_NVT::get_pressure() const
{
    // The pressure is given by
    // w / (3 * volume)
    // where w is
    // 48 * energy_scale * [r(i,j)^-12 - 0.5 * r(i,j)^-6] / n_particles
    // summed over distinct pairs (i,j), with r(i,j) being the distance
    // between particles i and j in reduced units.
    // In reduced units, the volume is 1 so the pressure becomes
    // 48 * w / 3.
    unsigned int n_coordinates = 3;
    std::vector<double> displacement(n_coordinates);
    double w(0), sq_distance;
    // Cycle over pairs of particles only.
    for(unsigned int i = 0; i < n_particles - 1; ++i)
        for(unsigned int j = i + 1; j < n_particles; ++j)
        {
            for(unsigned int d = 0; d < n_coordinates; ++d)
                displacement[d] = quotient(current_configuration[n_coordinates * i + d] - current_configuration[n_coordinates * j + d]);
            // In the computation of the potential en include only
            // the pairs whose distance is less than the cutoff radius.
            sq_distance = std::inner_product(displacement.begin(), displacement.end(), displacement.begin(), 0.);
            // In the expression for the potential energy the distance appears
            // in even powers only, so we don't need to compute the square
            // root.
            if(sq_distance < std::pow(distance_cutoff, 2))
                w += pow(sq_distance, -6) - 0.5 * pow(sq_distance, -3);
        }

    return 16. * w;
}

void metropolis_NVT::write_config(const std::string & output_file) const
{ 
    // Write the final configuration on file.
    std::cout << "Print final configuration to file " << output_file << "." << std::endl;
    std::ofstream output(output_file);

    // Output formatting.
    output.precision(4);
    output << std::scientific;
    const unsigned int col_width(16);

    unsigned int n_coordinates = 3;
    for(unsigned int i = 0; i < n_particles; ++i)
    {
        for(unsigned int d = 0; d < n_coordinates; ++d)
            // The lengths are output in units of the cell edge length.
            output << std::setw(col_width) << current_configuration[3 * i + d] / cell_edge_length;
        output << std::endl;
    }
    output.close();
    return;
}

void metropolis_NVT::write_config_xyz(const std::string & output_file) const
{ 
    // Write configuration in .xyz format.
    std::ofstream output(output_file);

    // Output formatting.
    output.precision(4);
    output << std::scientific;
    const unsigned int col_width(16);

    output << n_particles << "\n";
    output << "Simulation of a bunch of molecules interacting with a Lennard-Jones potential\n";
    unsigned int n_coordinates = 3;
    for(unsigned int i = 0; i < n_particles; ++i)
    {
        output << std::setw(col_width) << "LJ";
        for(unsigned int d = 0; d < n_coordinates; ++d)
            output << std::setw(col_width) << quotient(current_configuration[3 * i + d]);
        output << std::endl;
    }

    output.close();
}

double metropolis_NVT::quotient(double r) const
{  
    return r - cell_edge_length * std::round(r / cell_edge_length);
}

double metropolis_NVT::get_acceptance_rate() const
{
    return static_cast<double>(accepted_proposals) / total_proposals;
}
