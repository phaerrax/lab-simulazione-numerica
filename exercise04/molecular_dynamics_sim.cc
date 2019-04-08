#include "molecular_dynamics_sim.hh"
#include <sstream>
#include <numeric>
#include <iomanip>
#include <fstream>
#include <iostream>
#include <cmath>
// (string, vector and random.hh included in molecular_dynamics_sim.hh.)

molecular_dynamics_sim::molecular_dynamics_sim(const std::string & initial_config_file, double input_cell_edge_length)
{   
    cell_edge_length = input_cell_edge_length;
    n_particles = 0;

    // Clear the previous configuration (if there is any).
    position.clear();

    // Read the initial configuration from a file.
    std::cout << "Read initial configuration from file " << initial_config_file << "." << std::endl;
    std::ifstream input_config(initial_config_file);

    std::vector<double> new_point;
    std::string line, coordinate;

    // Read the file line by line.
    while(std::getline(input_config, line))
    {
        // Save each line in a new vector.
        // (The coordinates are given in units of the cell edge length, so
        // we need to rescale them.)
        new_point.clear();
        std::istringstream iss(line);
        while(iss >> coordinate)
            new_point.push_back(cell_edge_length * std::stod(coordinate));

        // Add the new vector to the list of positions of the particles.
        position.push_back(new_point);
        ++n_particles;
    }

    input_config.close();
    return;
}

molecular_dynamics_sim::molecular_dynamics_sim(const std::string & initial_config_file, const std::string & pre_initial_config_file, double input_cell_edge_length)
{
    cell_edge_length = input_cell_edge_length;
    n_particles = 0;

    // Clear the previous configuration (if there is any).
    position.clear();
    old_position.clear();

    // Read the initial configuration from a file.
    std::cout << "Read initial configuration from file " << initial_config_file << "." << std::endl;
    std::ifstream input_config(initial_config_file);

    std::vector<double> new_point;
    std::string line, coordinate;

    // Read the file line by line.
    while(std::getline(input_config, line))
    {
        // Save each line in a new vector.
        // (The coordinates are given in units of the cell edge length, so
        // we need to rescale them.)
        new_point.clear();
        std::istringstream iss(line);
        while(iss >> coordinate)
            new_point.push_back(cell_edge_length * std::stod(coordinate));

        // Add the new vector to the list of positions of the particles.
        position.push_back(new_point);
        ++n_particles;
    }

    input_config.close();

    // Same as above, but with the pre-initial configuration.
    std::cout << "Read pre-initial configuration from file " << pre_initial_config_file << "." << std::endl;
    input_config.open(pre_initial_config_file);

    // Read the file line by line.
    while(std::getline(input_config, line))
    {
        // Save each line in a new vector.
        // (The coordinates are given in units of the cell edge length, so
        // we need to rescale them.)
        new_point.clear();
        std::istringstream iss(line);
        while(iss >> coordinate)
            new_point.push_back(cell_edge_length * std::stod(coordinate));

        // Add the new vector to the list of positions of the particles.
        old_position.push_back(new_point);
    }

    input_config.close();

    // Calculate the initial velocity configuration.
    // 1. Compute the position at the first step using the Verlet algorithm,
    //    obtaining the current velocity too.
    // 2. Compute the temperature and compare it with the given one.
    // 3. Rescale the velocities to match the current temperature with the
    //    given one.
    // 4. Correct the old positions, using the newly computed velocities.

    // Before moving, save the original positions in a new vector.
    std::vector<std::vector<double>> original_init_position(position),
                                     original_preinit_position(old_position);

    move();
    // When the integration step is complete, the velocity of the particles
    // is calculated at the time instant corresponding to the initial
    // positions.

    // Restore the original position of the particles.
    position = original_init_position;
    old_position = original_preinit_position;

    return;
}

void molecular_dynamics_sim::set_temperature(double input_temperature)
{
    temperature = input_temperature;
    return;
}

void molecular_dynamics_sim::set_particle_number(unsigned int input_particle_number)
{
    n_particles = input_particle_number;
    return;
}

void molecular_dynamics_sim::set_particle_density(double input_particle_density)
{
    particle_density = input_particle_density;
    return;
}

void molecular_dynamics_sim::set_distance_cutoff(double input_distance_cutoff)
{
    distance_cutoff = input_distance_cutoff;
    return;
}

void molecular_dynamics_sim::set_integration_step(double input_integration_step)
{
    integration_step = input_integration_step;
    return;
}

void molecular_dynamics_sim::initialise_uniform(Random & rng)
{
    // Simulate an initial velocity configuration.
    std::cout << "Prepare random velocities with center of mass velocity equal to zero." << std::endl;

    unsigned int n_coordinates = position.begin()->size();
    // The number of dimensions, i.e. of coordinates of each velocity vector.

    // The velocities will be translated such that the velocity of the
    // centre of mass is zero. For this, we need the sum of the velocities.
    std::vector<double> com_velocity(n_coordinates, 0.);
    velocity.resize(n_particles);
    for(auto & v : velocity)
    {
        v.resize(n_coordinates);
        for(unsigned int d = 0; d < n_coordinates; ++d)
        {
            // Generate uniformly distributed velocities in [-0.5, 0.5).
            v[d] = rng.Rannyu() - 0.5;
            com_velocity[d] += v[d];
        }
    }

    for(unsigned int d = 0; d < n_coordinates; ++d)
        com_velocity[d] /= n_particles;

    // Subtract from each velocity the centre of mass velocity
    // calculated above.
    for(auto & v : velocity)
        for(unsigned int d = 0; d < n_coordinates; ++d)
            v[d] -= com_velocity[d];

    // From the initial position and the initial velocity, extrapolate
    // a value for the "pre-initial" position.
    old_position.resize(n_particles);
    for(unsigned int i = 0; i < n_particles; ++i)
    {
        old_position[i].resize(n_coordinates);
        for(unsigned int d = 0; d < n_coordinates; ++d)
            old_position[i][d] = position[i][d] - velocity[i][d] * time_step;
    }

    return;
}

void molecular_dynamics_sim::rescale_velocity(double input_temperature)
{
    temperature = input_temperature;
    // Compute the sum of the squares of the velocity of each particle.
    double sum_sq_velocities(0);
    for(auto & particle_velocity : velocity)
        sum_sq_velocities += std::inner_product(particle_velocity.begin(), particle_velocity.end(), particle_velocity.begin(), 0.);

    // Divide by the number of particles to obtain the mean square velocity.
    sum_sq_velocities /= n_particles;

    // Rescale all the velocities in order to match the input temperature
    // given before.
    unsigned int n_coordinates = position.begin()->size();
    double scale_factor = std::sqrt(3. * temperature / sum_sq_velocities);
    for(unsigned int i = 0; i < n_particles; ++i)
        for(unsigned int d = 0; d < n_coordinates; ++d)
        {
            velocity[i][d] *= scale_factor;
            // From the initial position and the initial velocity, recalculate
            // the value of the "pre-initial" position.
            old_position[i][d] = position[i][d] - velocity[i][d] * time_step;
        }

    return;
}

void molecular_dynamics_sim::initialise_maxwellboltzmann(double input_temperature, Random & rng)
{
    temperature = input_temperature;
    
    // Simulate an initial velocity configuration.
    std::cout << "Sample random velocities from a Maxwell-Boltzmann distribution." << std::endl;

    unsigned int n_coordinates = position.begin()->size();
    // The number of dimensions, i.e. of coordinates of each velocity vector.

    // Sampling...
    
    for(unsigned int i = 0; i < n_particles; ++i)
        for(unsigned int d = 0; d < n_coordinates; ++d)
            // From the initial position and the initial velocity, extrapolate
            // a value for the "pre-initial" position.
            old_position[i][d] = position[i][d] - velocity[i][d] * time_step;

    return;
}

void molecular_dynamics_sim::move(void)
{
    // Move particles using Verlet algorithm.
    std::vector<std::vector<double>> forces(n_particles);
    unsigned int n_coordinates = position.begin()->size();
    std::vector<double> new_position(n_coordinates);
    for(unsigned int i = 0; i < n_particles; ++i)
    { 
        for(unsigned int d = 0; d < n_coordinates; ++d)
        {
            // Knowing the current position, the previous position and the
            // force exerted on each particle, calculate the following point
            // in the particle trajectory.
            new_position[d] = 
                quotient(
                        2. * position[i][d] - old_position[i][d] + force(i, d) * std::pow(time_step, 2)
                        );
            // Calculate the current velocity.
            velocity[i][d] = quotient(new_position[d] - old_position[i][d]) / (2. * time_step);
        }

        // Move current position to old position, and new to current.
        old_position[i] = position[i];
        position[i] = new_position;
    }

    return;
}

double molecular_dynamics_sim::force(unsigned int this_particle, unsigned int dir) const
{
    // Compute forces as -grad V(r).
    double f(0), distance;
    unsigned int n_coordinates = position.begin()->size();
    std::vector<double> displacement(n_coordinates);

    for(unsigned int i = 0; i < n_particles; ++i)
    {
        if(i != this_particle)
        {
            for(unsigned int d = 0; d < n_coordinates; ++d)
                displacement[d] = quotient(position[this_particle][d] - position[i][d]);
            distance = std::sqrt(std::inner_product(displacement.begin(), displacement.end(), displacement.begin(), 0.));
            if(distance < distance_cutoff)
            {
                f += displacement[dir] * (48. * std::pow(distance, -14) - 24. * std::pow(distance, -8));
            }
        }
    }

    return f;
}

void molecular_dynamics_sim::measure() const
{
    // Compute thermodynamical quantities of the system, and append the
    // results to the given files.
    double potential_en(0),
           kinetic_en(0);

    std::ofstream potential_en_output,
                  kinetic_en_output,
                  temperature_output,
                  total_en_output;

    potential_en_output.open("output_epot.dat", std::ios::app);
    kinetic_en_output.open("output_ekin.dat", std::ios::app);
    temperature_output.open("output_temp.dat", std::ios::app);
    total_en_output.open("output_etot.dat", std::ios::app);

    double current_potential_en_density,
           current_kinetic_en_density,
           current_temperature,
           current_total_en_density;

    // Cycle over pairs of particles.
    unsigned int n_coordinates = position.begin()->size();
    std::vector<double> displacement(n_coordinates);
    double distance;
    for(unsigned int i = 0; i < n_particles - 1; ++i)
        for(unsigned int j = i + 1; j < n_particles; ++j)
        {
            for(unsigned int d = 0; d < n_coordinates; ++d)
                displacement[d] = quotient(position[i][d] - position[j][d]);
            // In the computation of the potential en include only
            // the pairs whose distance is less than the cutoff radius.
            distance = std::sqrt(std::inner_product(displacement.begin(), displacement.begin(), displacement.end(), 0.));
            if(distance < distance_cutoff)
                potential_en += 4. * pow(distance, -12) - 4. * pow(distance, -6);
        }          

    // Kinetic en (per the particle mass).
    for(auto & v : velocity)
        kinetic_en += 0.5 * std::inner_product(v.begin(), v.begin(), v.end(), 0.);

    current_potential_en_density = potential_en / n_particles;
    current_kinetic_en_density   = kinetic_en / n_particles;
    current_temperature          = 2. / 3. * kinetic_en / n_particles;
    current_total_en_density     = (kinetic_en + potential_en) / n_particles;

    potential_en_output << current_potential_en_density << std::endl;
    kinetic_en_output   << current_kinetic_en_density   << std::endl;
    temperature_output  << current_temperature          << std::endl;
    total_en_output     << current_total_en_density     << std::endl;

    potential_en_output.close();
    kinetic_en_output.close();
    temperature_output.close();
    total_en_output.close();

    return;
}

void molecular_dynamics_sim::write_config(const std::string & output_file) const
{ 
    // Write the final configuration on file.
    std::cout << "Print final configuration to file " << output_file << "." << std::endl;
    std::ofstream output(output_file);

    // Output formatting.
    output.precision(6);
    output << std::scientific;
    const unsigned int col_width(12);

    unsigned int n_coordinates = position.begin()->size();
    for(unsigned int i = 0; i < n_particles; ++i)
    {
        for(unsigned int d = 0; d < n_coordinates; ++d)
            // The lengths are output in units of the cell edge length.
            output << std::setw(col_width) << position[i][d] / cell_edge_length;
        output << std::endl;
    }
    output.close();
    return;
}

void molecular_dynamics_sim::write_config_xyz(const std::string & output_file_prefix, int n_conf) const
{ 
    // Write configuration in .xyz format.
    std::ofstream output(output_file_prefix + std::to_string(n_conf) + ".xyz");

    // The integer n_conf distinguishes between the "snapshots" of the system
    // taken at regular times during the simulation.
    // The list of config_N.xyz files will be read by an external program,
    // i.e. ovito, which will be able then to visualise the evolution of the
    // system.

    // Output formatting.
    output.precision(6);
    output << std::scientific;
    const unsigned int col_width(12);

    output << n_particles << "\n";
    output << "Simulation of a bunch of molecules interacting with a Lennard-Jones potential\n";
    unsigned int n_coordinates = position.begin()->size();
    for(unsigned int i = 0; i < n_particles; ++i)
    {
        output << std::setw(col_width) << "LJ";
        for(unsigned int d = 0; d < n_coordinates; ++d)
            output << std::setw(col_width) << quotient(position[i][d]);
        output << std::endl;
    }

    output.close();
}

double molecular_dynamics_sim::quotient(double r) const
{  
    return r - cell_edge_length * std::round(r / cell_edge_length);
}
