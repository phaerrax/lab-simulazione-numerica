#ifndef MOLECULAR_DYNAMICS_SIM
#define MOLECULAR_DYNAMICS_SIM

#include <string>
#include <vector>
#include <tuple>
#include "random.hh"

class molecular_dynamics_sim
{
    public:
        molecular_dynamics_sim(const std::string &, double = 1);
        molecular_dynamics_sim(const std::string &, const std::string &, double int_step, double = 1);
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

        // Extract thermodynamical quantities
        // ==================================
        double get_temperature();
        double get_potential_energy_density() const;
        double get_kinetic_energy_density();
        double get_pressure();
        // The temperature, pressure and kinetic energy methods, if possible,
        // use the same calculation of the mean square velocity, to save
        // some time.

		// Radial distribution of the particles
		// ====================================
		// Compute an histogram of the radial distribution of the particles of
		// the system.
		std::vector<double> get_radial_distribution(unsigned int) const;
		// Use the given number of bins.
		std::vector<double> get_radial_distribution(unsigned int, double) const;
		// Use the given number of bins, discarding the values over the given
		// threshold.

		// General measurement methods
		// ===========================
		std::tuple<double, double, double, double, std::vector<double>>
			measure(unsigned int, double);
		std::tuple<double, double, double, double, std::vector<double>>
			measure(unsigned int);
		// Measure everything at the same time. In order:
		// - temperature,
		// - potential energy,
		// - kinetic energy,
		// - pressure,
		// - radial distribution.
		// Calling the particular function for each quantity to be measured
		// requires looping over all (distinct) pairs of particles for every
		// function called, which wastes a lot of time.
		// Since the potential energy, the pressure and the radial
		// distribution all require a single loop, they can be effectively
		// calculated together.
		// The arguments are those relative to the radial distribution: the
		// number of bins and the maximum distance visualised in the
		// histogram (optional).

        // Output
        // ======
        void write_config(const std::string &) const;
        // Write the final configuration of the system in a file.
        void write_config_xyz(const std::string &) const;
        // Write the n-th configuration of the system in a file, in XYZ

        // format (to be read by ovito or similar programs).
        // This makes the visualisation of the progressive evolution of the
        // system as the algorithm advances possible.

    private:
		std::vector<double> force(unsigned int particle_index) const;
        // Calculate the force exerted on the particle at index
        // 'particle_index'.

        double quotient(double) const;
        std::vector<double> quotient(const std::vector<double> &) const;
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
               ms_velocity,
               total_volume,
               particle_density,
               cell_edge_length,
               integration_step,
               distance_cutoff;

        // Flags
        // =====
        bool ms_velocity_already_computed;
};

std::vector<double> operator-(const std::vector<double> &);
std::vector<double> operator+(const std::vector<double> &, const std::vector<double> &);
std::vector<double> operator-(const std::vector<double> &, const std::vector<double> &);

#endif
