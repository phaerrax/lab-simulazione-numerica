#include "metropolis_NVT.hh"
#include <algorithm>
#include <iterator>
#include <cmath>
#include <numeric>
#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>
#include <sstream>

metropolis_NVT::metropolis_NVT(const std::string & initial_config_file, double input_cell_edge_length):
    cell_edge_length(input_cell_edge_length),
    particle_density(0),
    accepted_proposals(0),
    total_proposals(0),
    current_configuration()
{   

    // Read the initial configuration from a file.
    std::cout << "Read initial configuration from file " << initial_config_file << "." << std::endl;
    std::ifstream input_config(initial_config_file);

    std::string line, coordinate;
    std::vector<double> new_point;

    // Read the file line by line, and save each line in a new vector.
    while(std::getline(input_config, line))
    {
        std::istringstream iss(line);
		// I expect that a line consists in all coordinates of a particle:
        // then the number of dimensions of the system is the number of
        // elements in a line.
		// The number of lines in the initial config file will be the
		// number of particles.
        new_point.clear();
        while(iss >> coordinate)
            // The coordinates are given in units of the cell edge length,
            // so I need to rescale them.
            new_point.push_back(cell_edge_length * std::stod(coordinate));

        current_configuration.push_back(new_point);
    }

    n_particles = current_configuration.size();
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

double metropolis_NVT::interaction_energy(unsigned int particle_n, const std::vector<double> & position) const
{
	double sq_distance,
	       energy(0);
    std::vector<double> displacement;
	for(unsigned int n = 0; n < n_particles; ++n)
		if(n != particle_n)
		{
            displacement = quotient(position - current_configuration[n]);
			sq_distance = std::inner_product(displacement.begin(), displacement.end(), displacement.begin(), 0.);
			if(sq_distance < std::pow(distance_cutoff, 2))
                // Lennard-Jones potential function.
				energy += std::pow(sq_distance, -6) - std::pow(sq_distance, -3);
		}

    return energy * 4;
}

void metropolis_NVT::next(Random & rng)
{
	// Select a particle at random.
	unsigned int selected_particle = static_cast<unsigned int>(rng.Rannyu(0, n_particles));
	std::vector<double> selected_position(current_configuration[selected_particle]);

	// Generate a proposal. 
	std::vector<double> proposed_position(selected_position);
    // Start from the previous position of the particle and move it with a
    // displacement sampled from a normal distribution.
    for(double & x : proposed_position)
		x = quotient(x + rng.Gauss(0, stdev));

    total_proposals++;

    // Calculate the energy of the current and proposed configurations.
	double energy_diff = interaction_energy(selected_particle, proposed_position) - 
        interaction_energy(selected_particle, selected_position);

	if(energy_diff < 0)
	{
		accepted_proposals++;
		current_configuration[selected_particle] = proposed_position;
	}
	else
	{
		double acceptance_threshold = std::exp(-energy_diff / temperature);
		if(rng.Rannyu() < acceptance_threshold)
		{
			accepted_proposals++;
			current_configuration[selected_particle] = proposed_position;
		}
	}
}

double metropolis_NVT::get_potential_energy_density() const
{
    std::vector<double> displacement;
    double potential_en(0), sq_distance;
    // Cycle over pairs of particles only.
    for(auto x = current_configuration.begin(); x != current_configuration.end(); ++x)
        for(auto y = x + 1; y != current_configuration.end(); ++y)
        {
            displacement = quotient(*x - *y);
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
    std::vector<double> displacement;
    double w(0), sq_distance;
    // Cycle over pairs of particles only.
    for(auto x = current_configuration.begin(); x != current_configuration.end(); ++x)
        for(auto y = x + 1; y != current_configuration.end(); ++y)
        {
            displacement = quotient(*x - *y);
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

    for(const auto & row : current_configuration)
    {
        for(auto x = row.begin(); x != row.end(); ++x)
            // The lengths are output in units of the cell edge length.
            output << std::setw(col_width) << *x / cell_edge_length;
        output << "\n";
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
    for(const auto & row : current_configuration)
    {
        output << std::setw(col_width) << "LJ";
        for(auto x = row.begin(); x != row.end(); ++x)
            output << std::setw(col_width) << quotient(*x);
        output << std::endl;
    }

    output.close();
}

double metropolis_NVT::quotient(double r) const
{  
    return r - cell_edge_length * std::round(r / cell_edge_length);
}

std::vector<double> metropolis_NVT::quotient(const std::vector<double> & r) const
{  
    std::vector<double> result;
    for(auto x : r)
        result.push_back(quotient(x));

    return result;
}

double metropolis_NVT::get_acceptance_rate() const
{
    return static_cast<double>(accepted_proposals) / total_proposals;
}

std::vector<double> operator+(const std::vector<double> & lhs, const std::vector<double> & rhs)
{
    std::vector<double> result;
    std::transform(lhs.begin(), lhs.end(), rhs.begin(), std::back_inserter(result), std::plus<double>());
    return result;
}

std::vector<double> operator-(const std::vector<double> & v)
{
    std::vector<double> result;
    std::transform(v.begin(), v.end(), std::back_inserter(result), std::negate<double>());
    return result;
}

std::vector<double> operator-(const std::vector<double> & lhs, const std::vector<double> & rhs)
{
    return lhs + (-rhs);
}
