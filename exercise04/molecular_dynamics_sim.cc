#include "molecular_dynamics_sim.hh"
#include <sstream>
#include <algorithm>
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
    ms_velocity_already_computed = false;

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

molecular_dynamics_sim::molecular_dynamics_sim(const std::string & initial_config_file, const std::string & pre_initial_config_file, double input_integration_step, double input_cell_edge_length):
    n_particles(0),
    cell_edge_length(input_cell_edge_length),
	integration_step(input_integration_step),
    ms_velocity_already_computed(false)
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

    // The move() method requires that the velocity vector is already of
    // the correct size.
    unsigned int n_coordinates = position.begin()->size();
    velocity.resize(n_particles);
    for(auto & v : velocity)
        v.resize(n_coordinates);

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
            old_position[i][d] = position[i][d] - velocity[i][d] * integration_step;
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
            old_position[i][d] = position[i][d] - velocity[i][d] * integration_step;
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

    // According to Boltzmann's theory, each component of the velocity
    // of a particle (under certain hypoteses) follows a normal distribution
    // with mean 0 and variance (2 * a)^-1, with
    // a = particle mass / (2 * Boltzmann constant * temperature).
    // In "Lennard-Jones units", the probability distribution is a normal
    // distribution with mean 0 and variance equal to the adimensional
    // temperature.
    velocity.resize(n_particles);
    std::vector<double> com_velocity(n_coordinates);
    for(auto & v : velocity)
    {
        v.resize(n_coordinates, 0);
        for(unsigned int d = 0; d < n_coordinates; ++d)
        {
            v[d] = rng.Gauss(0, std::sqrt(temperature));
            com_velocity[d] += v[d];
        }
    }

    // Correct the velocities so that there is no overall momentum.
    // Calculate the centre of mass velocity...
    for(unsigned int d = 0; d < n_coordinates; ++d)
        com_velocity[d] /= n_particles;

    // ...and subtract from each velocity the one of the centre of mass.
    for(auto & v : velocity)
        for(unsigned int d = 0; d < n_coordinates; ++d)
            v[d] -= com_velocity[d];
    
    old_position.resize(n_particles);
    for(unsigned int i = 0; i < n_particles; ++i)
    {
        old_position[i].resize(n_coordinates);
        for(unsigned int d = 0; d < n_coordinates; ++d)
            // From the initial position and the initial velocity, extrapolate
            // a value for the "pre-initial" position.
            old_position[i][d] = position[i][d] - velocity[i][d] * integration_step;
    }

    return;
}

void molecular_dynamics_sim::move(void)
{
    ms_velocity_already_computed = false;
    // As the algorithm moves to the next step, the previous value of the
    // mean square velocity is no longer valid since it refers to the
    // a past state of the system.

    // Move particles using Verlet algorithm.
    unsigned int n_coordinates = position.begin()->size();
    std::vector<double> new_position(n_coordinates);

	/*
	   Be careful not to merge this loop with the other one below: they must
	   be separated because the forces are supposed to be calculated using the
	   current positions of the particles.
	   Otherwise, the 2nd particle will have its new position calculated using
	   the already updated position of the 1st particle, and so on.
	*/
	std::vector<std::vector<double>> f(n_particles);
    for(unsigned int i = 0; i < n_particles; ++i)
		f[i] = force(i);

    for(unsigned int i = 0; i < n_particles; ++i)
    { 
        for(unsigned int d = 0; d < n_coordinates; ++d)
        {
            // Knowing the current position, the previous position and the
            // force exerted on each particle, calculate the following point
            // in the particle trajectory:
            // x(t + s) = 2 * x(t) - x(t - s) + s^2 * f(x(t)).
            new_position[d] = 
                quotient(
                        2. * position[i][d] - old_position[i][d] + f[i][d] * std::pow(integration_step, 2)
                        );
            // Calculate the current velocity:
            // v(t) = (x(t + s) - x(t - s)) / 2.
            velocity[i][d] = quotient(new_position[d] - old_position[i][d]) / (2. * integration_step);
        }

        // Move current position to old position, and new to current.
        old_position[i] = position[i];
        position[i] = new_position;
    }

    return;
}

std::vector<double> molecular_dynamics_sim::force(unsigned int n) const
{
	// 0 <= n < n_particles.
    // Compute (adimensional) forces as -grad V(x).
    double sq_distance;
	auto this_particle = position.begin() + n;
    std::vector<double> displacement,
						f(position.begin()->size());
    for(auto x = position.begin(); x != position.end(); ++x)
    {
        if(x != this_particle)
        {
			displacement = quotient(*this_particle - *x);
            sq_distance = std::inner_product(displacement.begin(), displacement.end(), displacement.begin(), 0.);
            // In the expression for the force the distance appears in even
            // powers only, so we don't need to compute the square root.
            if(sq_distance < std::pow(distance_cutoff, 2))
			{
				for(unsigned int d = 0; d < displacement.size(); ++d)
					f[d] += displacement[d] * 24 * (2 * std::pow(sq_distance, -7) - std::pow(sq_distance, -4));
			}
        }
    }

    return f;
}

double molecular_dynamics_sim::get_temperature() const
{
    if(!ms_velocity_already_computed)
    {
        ms_velocity = 0;
        for(auto & v : velocity)
            ms_velocity += std::inner_product(v.begin(), v.end(), v.begin(), 0.);
        ms_velocity_already_computed = true;
    }

    return ms_velocity / (3 * n_particles);
}

double molecular_dynamics_sim::get_potential_energy_density() const
{
    std::vector<double> displacement;
    double potential_en(0), sq_distance;
    // Cycle over pairs of particles only.
    for(auto x = position.begin(); x != position.end(); ++x)
        for(auto y = x + 1; y != position.end(); ++y)
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

double molecular_dynamics_sim::get_kinetic_energy_density() const
{
    if(!ms_velocity_already_computed)
    {
        ms_velocity = 0;
        for(auto & v : velocity)
            ms_velocity += std::inner_product(v.begin(), v.end(), v.begin(), 0.);
        ms_velocity_already_computed = true;
    }

    return ms_velocity / (2 * n_particles);
}

double molecular_dynamics_sim::get_pressure() const
{
    // We can save some time if the temperature has already been calculated.
    if(!ms_velocity_already_computed)
    {
		/*
		   This is called ms_velocity but really it is the sum
		   of the square velocity of each particle. There is
		   nothing "mean" here. I should probably change the
		   name of this variable.
		   The code treats this variable correctly, though, so
		   there should be no bugs related to this.
		*/
        ms_velocity = 0;
        for(auto & v : velocity)
            ms_velocity += std::inner_product(v.begin(), v.end(), v.begin(), 0.);
        ms_velocity_already_computed = true;
    }

    // The pressure is given by
    // density * Boltzmann constant * temperature + w / (3 * volume)
    // where w is
    // 48 * energy_scale * [r(i,j)^-12 - 0.5 * r(i,j)^-6] / n_particles
    // summed over distinct pairs (i,j), with r(i,j) being the distance
    // between particles i and j in reduced units.
    // In reduced units, the volume is 1 so the pressure becomes
    // density * temperature + w / 3.
    // (In the calculations below we factored out the 48.)
    std::vector<double> displacement;
    double w(0), sq_distance;
    // Cycle over pairs of particles only.
    for(auto x = position.begin(); x != position.end(); ++x)
		for(auto y = x + 1; y != position.end(); ++y)
		{
			displacement = quotient(*x - *y);
            sq_distance = std::inner_product(displacement.begin(), displacement.end(), displacement.begin(), 0.);
            // In the computation of the pressure include only
            // the pairs whose distance is less than the cutoff radius.
            if(sq_distance < std::pow(distance_cutoff, 2))
                w += std::pow(sq_distance, -3) * (std::pow(sq_distance, -3) - 0.5);
        }

    return particle_density * ms_velocity / (3. * n_particles) + 16. * w * std::pow(cell_edge_length, -3);
}

std::tuple<double, double, double, double>
molecular_dynamics_sim::measure() const
{
    std::vector<double> displacement;

    double temperature,
		   distance,
		   potential_en(0),
		   kinetic_en,
		   pressure,
		   w(0);

    // Cycle over pairs of particles only.
    for(auto x = position.begin(); x != position.end(); ++x)
        for(auto y = x + 1; y != position.end(); ++y)
        {
            displacement = quotient(*x - *y);
            distance = std::sqrt(std::inner_product(displacement.begin(), displacement.end(), displacement.begin(), 0.));

			// Potential energy and pressure.
            if(distance < distance_cutoff)
			{
                potential_en += 4. * std::pow(distance, -12) - 4. * std::pow(distance, -6);
                w += std::pow(distance, -12) - 0.5 * std::pow(distance, -6);
			}

		}

    potential_en /= n_particles;

	ms_velocity = 0;
	for(auto & v : velocity)
		ms_velocity += std::inner_product(v.begin(), v.end(), v.begin(), 0.);
	ms_velocity_already_computed = true;

    temperature = ms_velocity / (3. * n_particles);
    kinetic_en = ms_velocity / (2 * n_particles);
    pressure = particle_density * temperature + 16. * w * std::pow(cell_edge_length, -3);

	return std::make_tuple(
			temperature,
			potential_en,
			kinetic_en,
			pressure
			);
}

void molecular_dynamics_sim::write_config(const std::string & output_file) const
{ 
    // Write the final configuration on file.
    std::cout << "Print final configuration to file " << output_file << "." << std::endl;
    std::ofstream output(output_file);

    // Output formatting.
    output.precision(4);
    output << std::scientific;
    const unsigned int col_width(16);

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

void molecular_dynamics_sim::write_config_xyz(const std::string & output_file) const
{ 
    // Write configuration in .xyz format.
    std::ofstream output(output_file);

    // Output formatting.
    output.precision(4);
    output << std::scientific;
    const unsigned int col_width(16);

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
    return r - cell_edge_length * rint(r / cell_edge_length);
}

std::vector<double> molecular_dynamics_sim::quotient(const std::vector<double> & r) const
{  
    std::vector<double> result;
    for(auto x : r)
        result.push_back(quotient(x));

    return result;
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
