#ifndef METROPOLIS_NVT_HH
#define METROPOLIS_NVT_HH

#include "random.hh"
#include <vector>
#include <string>

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
        double total_volume,
               particle_density,
			   temperature,
               cell_edge_length,
               distance_cutoff,
			   current_potential_energy,
			   stdev;

        unsigned int n_particles; // Total number of particles in the system.
		unsigned int accepted_proposals, total_proposals;
        
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
