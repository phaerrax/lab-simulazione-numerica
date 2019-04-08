#ifndef MOLECULAR_DYNAMICS_SIM
#define MOLECULAR_DYNAMICS_SIM

#include <string>
#include <vector>
#include "random.hh"

class molecular_dynamics_sim
{
    public:
        molecular_dynamics_sim(const std::string &, double = 1);
        molecular_dynamics_sim(const std::string &, const std::string &, double = 1);
        // Construct the class, collecting from a file the initial
        // configuration of the system. The length of the cell edge can be
        // supplied in the last argument, in order to rescale the coordinates
        // from the file.
        // An additional file can be provided, containing the positions of
        // the particles in a time step previous to the initial instant.
        
        // If the pre-initial configuration is provided, the constructor also
        // calculates the position at the next step with the Verlet algorithm,
        // then in computes the velocity of the particles.

        molecular_dynamics_sim() = delete;
        // Simulating the evolution of the system without an initial
        // configuration is meaningless.

        // Set internal parameters
        // =======================
        void set_temperature(double);
        void set_particle_number(unsigned int);
        void set_particle_density(double);
        void set_distance_cutoff(double);
        void set_integration_step(double);

        // Velocity initialisation
        // =======================
        void initialise_uniform(Random & rng);
        // Generate an initial distribution of velocities, from which
        // simulate the position of the particles in the time step just
        // before the start.
        void initialise_maxwellboltzmann(double, Random & rng);
        // Same as above, except this time the velocities are sampled from
        // a Maxwell-Boltzmann distribution in accord with the given
        // temperature. The temperature is set to the given value.

        void rescale_velocity(double);
        // If the velocities of the particles have been generated with methods
        // such as initialise_uniform() or the two-file constructor, the
        // temperature will depend the velocity distribution.
        // One may want to set the system to a desired temperature: this
        // method rescales the velocity of the particles in order to match
        // the desired temperature. This also recalculates the old positions
        // from the current position and the newly computed velocities,
        // and sets the temperature to the given value.

        // Knowing the position at the current and previous steps, and the
        // force exerted on each particle, the algorithm can start.

        void move();
        // Advance the algorithm one step with the Verlet method.

        double measure() const;
        // Read some thermodynamic properties of the system, calculated from
        // its current state.
        // In order not to use too many resources, the thermodynamical
        // properties of the system should not be evaluated at each new step,
        // but only after a certain number of steps, say every tenth step.
        // This has to be done directly in the program using the class,
        // by timing the calls to move() and measure().

        // Output
        // ======
        void write_config(const std::string &) const;
        // Write the final configuration of the system in a file.
        void write_config_xyz(const std::string &, int) const;
        // Write the n-th configuration of the system in a file, in XYZ
        // format (to be read by ovito or similar programs).
        // This makes the visualisation of the progressive evolution of the
        // system as the algorithm advances possible.

    private:
        double force(unsigned int particle_index, unsigned int dir) const;
        // Calculate the force exerted on the particle at index
        // 'particle_index', in the direction 'dir'.

        double quotient(double) const;
        // Quotient the position/distance/etc. in order to reduce everything
        // to the unit cell of the lattice.
        // This imposes periodic boundary conditions on the unit cell of
        // the lattice.

        std::vector<std::vector<double>> position, old_position, velocity;
        // Position and veolcity of the particles: each sub-array is a
        // d-tuple containing the position/velocity of the n-th particle.
        // Another array contains the position of the particles at the
        // previous integration step.

        // Internal / initial parameters
        unsigned int n_particles; // Total number of particles in the system.
        double temperature,
               total_volume,
               time_step,
               particle_density,
               cell_edge_length,
               integration_step,
               distance_cutoff;
};

#endif
