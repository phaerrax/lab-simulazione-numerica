#include "molecular_dynamics_sim.hh"
#include <string>
#include <sstream>
#include <numeric>
#include <iomanip>
#include "random.hh"

void molecular_dynamics_sim::set_init_config(const std::string & input_parameters_file)
{   
    // It is assumed that the random number generator is already initialised.

    double ep, ek, pr, et, vir; // ??

    std::cout << "Classic Lennard-Jones fluid\n
                  Molecular dynamics simulation in NVE ensemble\n\n
                  Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]\n\n
                  The program uses Lennard-Jones units" << std::endl;

    // Read the input parameters from a file.
    std::cout << "Read input parameters from file " << input_parameters_file << "." << std::endl;
    std::ifstream input_parameters(input_parameters_file);

    input_parameters >> input_temperature;

    input_parameters >> n_particles;
    std::cout << "Number of particles: " << n_particles << "\n";

    input_parameters >> particle_density;
    std::cout << "Density of particles: " << particle_density << "\n";

    total_volume = static_cast<double>(n_particles) / particle_density;
    std::cout << "Volume of the simulation cell: " << total_volume << "\n";

    cell_edge_length = std::pow(total_volume, 1. / 3.);
    std::cout << "Edge of the simulation cell: " << cell_edge_length << "\n";

    input_parameters >> distance_cutoff;
    input_parameters >> time_step;
    input_parameters >> n_steps;
    input_parameters >> print_steps;

    std::cout << "The program integrates Newton equations with the Verlet method, using\n
                  Time step: " << time_step << ",\n
                  Number of steps: " << n_steps << "." << std::endl;

    input_parameters.close();
    return;
}

void molecular_dynamics_sim::initialise_uniform(const std::string & initial_config_file, Random & rng)
{
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
        std::istringstream iss(line);
        while(iss >> coordinate)
            new_point.push_back(cell_edge_length * std::stod(coordinate));

        // Add the new vector to the list of positions of the particles.
        position.push_back(new_point);
    }
    // Maybe check that the length of the 'position' vector actually equals
    // the variable 'n_particles'.

    input_config.close();

    // Simulate an initial velocity configuration.
    std::cout << "Prepare random velocities with center of mass velocity equal to zero." << std::endl;

    unsigned int n_coordinates = position.begin()->size();
    // The number of dimensions, i.e. of coordinates of each velocity vector.

    std::vector<double> sum(n_coordinates, 0.);
    for(auto & particle_velocity : velocity)
        for(unsigned int d = 0; d < n_coordinates; ++d)
        {
            // Generate uniformly distributed velocities in [-0.5, 0.5).
            particle_velocity[d] = rng.Rannyu() - 0.5;
            // The velocities will be translated such that the velocity of
            // the centre of mass is zero. For this, we need the sum of
            // the velocities.
            sum[d] += particle_velocity[d];
        }

    for(unsigned int i = 0; i < n_coordinates; ++i)
        sum[i] /= n_particles;

    double sum_sq_velocities;
    for(auto & particle_velocity : velocity)
        for(unsigned int i = 0; i < n_coordinates; ++i)
        {
            // Subtract from each velocity the centre of mass velocity
            // calculated above.
            particle_velocity[i] -= sum[i];
            // Compute the sum of the squares of the velocity of each
            // particle.
            sum_sq_velocities += std::inner_product(particle_velocity.begin(), particle_velocity.end(), particle_velocity.begin(), 0.);
        }

    // Divide by the number of particles to obtain the mean square velocity.
    sum_sq_velocities /= n_particles;

    // Rescale all the velocities in order to match the input temperature
    // given before.
    double scale_factor = std::sqrt(3. * input_temperature / sum_sq_velocities);
    for(auto & particle_velocity : velocity)
        for(unsigned int i = 0; i < n_coordinates; ++i)
        {
            particle_velocity[i] *= scale_factor;
            // From the initial position and the initial velocity, extrapolate
            // a value for the "pre-initial" position.
            old_position[i] = position[i] - particle_velocity[i] * time_step;
    }

    return;
}

void molecular_dynamics_sim::initialise_maxwellboltzmann(const std::string & initial_config_file, Random & rng)
{
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
    }
    // Maybe check that the length of the 'position' vector actually equals
    // the variable 'n_particles'.

    input_config.close();

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

void molecular_dynamics_sim::initialise_from_file(const std::string & initial_config_file, const std::string & pre_initial_config_file)
{
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
    }
    // Maybe check that the length of the 'position' vector actually equals
    // the variable 'n_particles'.

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
    std::vector<std::vector<double>> original_init_position(position);

    move();
    double sum_sq_velocities;
    for(auto & particle_velocity : velocity)
        // Compute the sum of the squares of the velocity of each particle.
        sum_sq_velocities += std::inner_product(particle_velocity.begin(), particle_velocity.end(), particle_velocity.begin(), 0.);

    // Divide by the number of particles to obtain the mean square velocity.
    sum_sq_velocities /= n_particles;

    // Rescale all the velocities in order to match the input temperature
    // given before.
    double scale_factor = std::sqrt(3. * input_temperature / sum_sq_velocities);
    for(unsigned int i = 0; i < n_particles; ++i) 
        for(unsigned int d = 0; d < n_coordinates; ++d)
        {
            velocity[i][d] *= scale_factor;

            // From the (original) initial position and the initial
            // velocity, extrapolate a value for the "pre-initial" position.
            old_position[i][d] = original_init_position[i][d] - velocity[i][d] * time_step;
        }
    // Restore the original values (as read from the file) of the initial
    // positions.
    position = original_init_position;

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
            // Verlet integration scheme.
            new_position[d] = 
                quotient(
                        2. * position[i][d] - old_position[i][d] + force(i, d) * std::pow(time_step, 2)
                        );
            velocity[i][d] = quotient(new_position[d] - old_position[i][d]) / (2. * time_step);
        }

        // Move current position to old position, and new to current.
        old_position[i] = position[i];
        position[i] = new_position;
    }
    return;
}

double molecular_dynamics_sim::force(unsigned int this_particle, unsigned int dir)
{
    // Compute forces as -grad V(r)
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
            if(dr < distance_cutoff)
            {
                f += displacement[d] * (48. * std::pow(distance, -14) - 24. * std::pow(distance, -8));
            }
        }
    }

    return f;
}

void molecular_dynamics_sim::measure()
{
    // Compute thermodynamical quantities of the system, and append the
    // results to the given files.
    double potential_en(0),
           kinetic_en(0),
           pair_interaction_en;

    std::ofstream potential_en_output,
                  kinetic_en_output,
                  temperature_output,
                  total_en_output;

    potential_en_output("output_epot.dat", std::ios::app);
    kinetic_en_output("output_ekin.dat", std::ios::app);
    temperature_output("output_temp.dat", std::ios::app);
    total_en_output("output_etot.dat", std::ios::app);

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
    for(auto & v : velocities)
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
    std::ofstream output_file(output_file);

    // Output formatting.
    output_file.precision(6);
    output_file << std::scientific;
    const unsigned int col_width(12);

    unsigned int n_coordinates = position.begin()->size();
    for(unsigned int i = 0; i < n_particles; ++i)
    {
        for(unsigned int d = 0; d < n_coordinates; ++d)
            // The lengths are output in units of the cell edge length.
            output_file << std::setw(col_width) << position[i][d] / cell_edge_length;
        output_file << std::endl;
    }
    output_file.close();
    return;
}

void molecular_dynamics_sim::write_config_xyz(const std::string & output_file, int n_conf) const
{ 
    // Write configuration in .xyz format.
    std::ofstream output_file("frames/config_" + std::to_string(nconf) + ".xyz");

    // The integer n_conf distinguishes between the "snapshots" of the system
    // taken at regular times during the simulation.
    // The list of config_N.xyz files will be read by an external program,
    // i.e. ovito, which will be able then to visualise the evolution of the
    // system.

    // Output formatting.
    output_file.precision(6);
    output_file << std::scientific;
    const unsigned int col_width(12);

    output_file << n_particles << "\n";
    output_file << "Simulation of a bunch of molecules interacting with a Lennard-Jones potential\n";
    for(unsigned int i = 0; i < n_particles; ++i)
    {
        output_file << "LJ  " << quotient(x[i]) << "   " <<  quotient(y[i]) << "   " << quotient(z[i]) << std::endl;
    }
    unsigned int n_coordinates = position.begin()->size();
    for(unsigned int i = 0; i < n_particles; ++i)
    {
        output_file << std::setw(col_width) << "LJ";
        for(unsigned int d = 0; d < n_coordinates; ++d)
        {
            output_file << std::setw(col_width) << quotient(position[i][d])
        }
        output_file << std::endl;
    }

    output_file.close();
}

double molecular_dynamics_sim::quotient(double r) const
{  
    return r - cell_edge_length * std::round(r / cell_edge_length);
}
