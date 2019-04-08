#ifndef MOLECULAR_DYNAMICS_SIM
#define MOLECULAR_DYNAMICS_SIM

class molecular_dynamics_sim
{
    public:
        void set_init_config(const std::string &);
        // Collect from a file the initial configuration of the system.

        void initialise_uniform(const std::string &, Random & rng);
        // Generate an initial distribution of velocities, from which
        // simulate the position of the particles in the time step just
        // before the start.
        // The velocities are drawn from a uniform distribution, then
        // they are rescaled in order to match them with a given temperature.
        // Knowing the position at the current and previous steps, and the
        // acceleration exerted on each particle, the algorithm can start.
        void initialise_maxwellboltzmann(const std::string &, Random & rng);
        // Same as above, except this time the velocities are sampled from
        // a Maxwell-Boltzmann distribution in accord with the given
        // temperature (this eliminates the need to rescale the velocities
        // afterwards).
        void initialise_from_file(const std::string &, const std::string &);
        // The "pre-initial" configuration can be also be supplied as a
        // text file; this way the position at the next step can be
        // calculated directly with the Verlet algorithm.
        // Then the velocity of the particles can be calculated, and it is
        // rescaled in order to match the previously set temperature.

        void move();
        // Advance the algorithm one step with the Verlet method.

        // Read some thermodynamic properties of the system, calculated from
        // its current state.
        double measure() const;

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
        // triplet containing the position/velocity of the n-th particle.
        // Another array contains the position of the particles at the
        // previous integration step.

        // Internal / initial parameters
        unsigned int n_particles, // Total number of particles in the system.
                     n_steps,     // Number of integration steps.
                     print_steps;
        // In order not to use too many resources, the thermodynamical
        // properties of the system are not evaluated at each new step, but
        // only after 'print_steps' steps.
        double input_temperature,
               total_volume,
               time_step,
               particle_density,
               cell_edge_length,
               distance_cutoff;

        // ??
        double potential_energy,
               kinetic_energy,
               total_energy,
               temperature;

};

#endif
