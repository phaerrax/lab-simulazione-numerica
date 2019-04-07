#include "molecular_dynamics_sim.hh"
#include <string>
#include <sstream>
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
        for(unsigned int i = 0; i < n_coordinates; ++i)
        {
            // Generate uniformly distributed velocities in [-0.5, 0.5).
            particle_velocity[i] = rng.Rannyu() - 0.5;
            // The velocities will be translated such that the velocity of
            // the centre of mass is zero. For this, we need the sum of
            // the velocities.
            sum[i] += particle_velocity[i];
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
    double potential_energy(0),
           kinetic_energy(0),
           vij; // ??
    double dx, dy, dz, dr;
    std::ofstream Epot, Ekin, Etot, Temp;

    Epot.open("output_epot.dat", std::ios::app);
    Ekin.open("output_ekin.dat", std::ios::app);
    Temp.open("output_temp.dat", std::ios::app);
    Etot.open("output_etot.dat", std::ios::app);

    v = 0.0; //reset observables
    t = 0.0;

    //cycle over pairs of particles
    for(unsigned int i = 0; i < n_particles-1; ++i)
    {
        for(int j=i+1; j<n_particles; ++j)
        {
            dx = quotient( x[i] - x[j] );
            dy = quotient( y[i] - y[j] );
            dz = quotient( z[i] - z[j] );

            dr = dx*dx + dy*dy + dz*dz;
            dr = sqrt(dr);
            if(dr < distance_cutoff)
            {
                vij = 4.0/pow(dr,12) - 4.0/pow(dr,6);
                //Potential energy
                v += vij;
            }
        }          
    }

    //Kinetic energy
    for(unsigned int i = 0; i < n_particles; ++i)
        t += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);

    stima_pot = v/static_cast<double>(n_particles); //Potential energy
    stima_kin = t/static_cast<double>(n_particles); //Kinetic energy
    stima_temp = (2.0 / 3.0) * t/static_cast<double>(n_particles); //Temperature
    stima_etot = (t+v)/static_cast<double>(n_particles); //Total enery

    Epot << stima_pot  << std::endl;
    Ekin << stima_kin  << std::endl;
    Temp << stima_temp << std::endl;
    Etot << stima_etot << std::endl;

    Epot.close();
    Ekin.close();
    Temp.close();
    Etot.close();

    return;
}

void molecular_dynamics_sim::write_config(const std::string & output_file) const
{ 
    // Write the final configuration on file.
    std::cout << "Print final configuration to file " << output_file << "." << std::endl;
    std::ofstream output_file(output_file);

    for(unsigned int i = 0; i < n_particles; ++i)
    {
        output_file << x[i]/cell_edge_length << "   " <<  y[i]/cell_edge_length << "   " << z[i]/cell_edge_length << std::endl;
    }
    output_file.close();
    return;
}

void molecular_dynamics_sim::write_config_xyz(const std::string & output_file, int nconf) const
{ 
    // Write configuration in .xyz format.
    std::ofstream output_file;

    output_file.open("frames/config_" + std::to_string(nconf) + ".xyz");
    output_file << n_particles << std::endl;
    output_file << "This is only a comment!" << std::endl;
    for(unsigned int i = 0; i < n_particles; ++i)
    {
        output_file << "LJ  " << quotient(x[i]) << "   " <<  quotient(y[i]) << "   " << quotient(z[i]) << std::endl;
    }
    output_file.close();
}

double molecular_dynamics_sim::quotient(double r) const
{  
    return r - cell_edge_length * std::round(r / cell_edge_length);
}
