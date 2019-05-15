#ifndef METROPOLIS_NVT_HH
#define METROPOLIS_NVT_HH

#include "random.hh"
#include <vector>
#include <string>
#include <tuple>

class metropolis_NVT
{
    public:
        metropolis_NVT(const std::string &, double = 1);
        // Construct the class, collecting from a file the initial
        // configuration of the system. The length of the cell edge can be
        // supplied in the last argument, in order to rescale the coordinates
        // from the file.

        metropolis_NVT() = delete;
        // Simulating the evolution of the system without an initial
        // configuration is meaningless.
        
        // Generate a new configuration, then accept it or keep
        // the current one.
        void next(Random &);

        // Set internal parameters
        // =======================
        void set_particle_density(double);
        void set_temperature(double);
        void set_distance_cutoff(double);
        void set_step_stdev(double);

        // Extract thermodynamical quantities
        // ==================================
        double get_potential_energy_density() const;
        double get_pressure() const;
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
		std::tuple<double, double, std::vector<double>>
			measure(unsigned int, double) const;
		std::tuple<double, double, std::vector<double>>
			measure(unsigned int) const;
		// Measure everything at the same time. In order:
		// - potential energy,
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

        // Utilities
        double size() const;
        double get_acceptance_rate() const;

    private:
        double interaction_energy(unsigned int, const std::vector<double> &) const;
        // Calculates the interaction energy of the particle at the given index, if it were at the given positon.

        double quotient(double) const;
        std::vector<double> quotient(const std::vector<double> &) const;
        // Quotient the position/distance/etc. in order to reduce everything
        // to the unit cell of the lattice.
        // This imposes periodic boundary conditions on the unit cell of
        // the lattice.

        // Internal parameters
        double cell_edge_length,
               particle_density,
			   temperature,
               distance_cutoff,
			   current_potential_energy,
			   stdev;

        unsigned int n_particles, // Total number of particles in the system.
		             accepted_proposals, total_proposals;
        
        // Holds the current configuration. 
        std::vector<std::vector<double>> current_configuration;
        // I could use arrays instead of vectors, and define this class as a
        // template on the number of points, since I believe there are
        // no use cases in which the number of points should change during the
        // simulation, but:
        // - the number of points may very well be specified in a sort of
        //   "input file" therefore it may not be known at compile-time;
        // - the number of points can be very large, and other arrays of the
        //   same size would be created during the simulation steps, occupying
        //   a lot of space in the stack (this may be an exaggeration...).
};

std::vector<double> operator-(const std::vector<double> &);
std::vector<double> operator+(const std::vector<double> &, const std::vector<double> &);
std::vector<double> operator-(const std::vector<double> &, const std::vector<double> &);

#endif
