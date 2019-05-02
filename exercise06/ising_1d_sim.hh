#ifndef ISING_1D_SIM
#define ISING_1D_SIM

#include "random.hh"
#include <vector>

class ising_1d_sim
{
	public:
		ising_1d_sim(unsigned int, Random & rng);
		// Defines a one-dimensional Ising model on a circle, i.e. a system
		// of the given number of sites with periodic boundary conditions.
		// The sites will be filled with the random number generator suppled
		// as an argument.
	
        // Set internal parameters
        // =======================
        void set_inverse_temperature(double);
        void set_interaction_energy(double);
        void set_external_magnetic_field(double);

		// Advance to another state
		// ========================
		void next_metropolis(Random &);
		void next_gibbs(Random &);

		// Extract physical quantities
		// ===========================
		double get_energy() const;
		// Returns the value of the Hamiltonian for the current spin
		// configuration.
		int get_spin_sum() const;
		// Returns the sum of all spins of the current configuration.

		// Get internal parameters
		// =======================
		unsigned int size() const;
	private:
		std::vector<int> spins;

		double inverse_temperature,
			   interaction_energy,
			   ext_magnetic_field;

        int quotient(int) const;
        // Impose periodic boundary conditions on the unit cell of
        // the lattice.
};

#endif
