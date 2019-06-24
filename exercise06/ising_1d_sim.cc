#include "ising_1d_sim.hh"
#include "random.hh"
#include <numeric>
#include <numeric>
#include <algorithm>
#include <cmath>

ising_1d_sim::ising_1d_sim(unsigned int n_spins, Random & rng)
{
	spins.resize(n_spins);
	inverse_temperature = 0;
	interaction_energy = 0;
	ext_magnetic_field = 0;
	
	for(auto & s : spins)
	{
		if(rng.Rannyu() < 0.5)
			s = -1;
		else
			s = 1;
	}
	return;
}

void ising_1d_sim::set_inverse_temperature(double input)
{
	inverse_temperature = input;
	return;
}

void ising_1d_sim::set_interaction_energy(double input)
{
	interaction_energy = input;
	return;
}

void ising_1d_sim::set_external_magnetic_field(double input)
{
	ext_magnetic_field = input;
	return;
}

void ising_1d_sim::next_metropolis(Random & rng)
{
	// Try to flip every spin in the configuration with the Metropolis
	// algorithm.
	double energy_diff;
	std::vector<unsigned int> draw_order(spins.size());
	std::iota(draw_order.begin(), draw_order.end(), 0);
	// Select, in random order, all the spins of the system.
	std::random_shuffle(
			std::begin(draw_order),
			std::end(draw_order),
			// std::random_shuffle needs a function that generates (uniformly)
			// a non-negative value less than its argument.
			[&rng](int n) {
				return static_cast<int>(rng.Rannyu(0, n));
			}
			);

	for(auto i : draw_order)
	{
		energy_diff = 2 * spins[quotient(i)] * (interaction_energy * (spins[quotient(i + 1)] + spins[quotient(i - 1)]) + ext_magnetic_field);
        // The acceptance threshold is less than 1 iff the energy difference
        // is greater than zero; this means that if the energy difference is
        // less than zero the step is automatically accepted.
		if(energy_diff > 0.)
		{
			if(rng.Rannyu() < std::exp(-inverse_temperature * energy_diff))
				spins[i] *= -1;
		}
		else
			spins[i] *= -1;
	}

	return;
}

void ising_1d_sim::next_gibbs(Random & rng)
{
	// Try to flip all the spins, sequentially, using the probability of the
	// new spin value conditioned on all the other values.
	double x;
	for(unsigned int i = 0; i < spins.size(); ++i)
	{
		x = inverse_temperature * (interaction_energy * (spins[quotient(i + 1)] + spins[quotient(i - 1)]) + ext_magnetic_field);
		// The probability that the spin will be -1 is 1 / (1 + e^(2 * x)).
		if(rng.Rannyu() * (1. + std::exp(2. * x)) < 1.)
			spins[i] = -1;
		else
			spins[i] = 1;
	}
	return;
}

double ising_1d_sim::get_energy() const
{
	double H(0);
	for(unsigned int i = 0; i < spins.size(); ++i)
		H += -interaction_energy * spins[i] * spins[quotient(i + 1)] - 0.5 * ext_magnetic_field * (spins[i] + spins[quotient(i + 1)]);

	return H;
}

int ising_1d_sim::get_spin_sum() const
{
	return std::accumulate(spins.begin(), spins.end(), 0);
}

unsigned int ising_1d_sim::size() const
{
	return spins.size();
}

int ising_1d_sim::quotient(int i) const
{
	// The sites are numbered from 0 to n_spins - 1; if the index is out of these
	// bounds, reduce it to its equivalent index inside the bounds (mod n_spins).
	int n_spins(spins.size());
    if(i >= n_spins)
		i -= n_spins;
    else
		if(i < 0)
			i += n_spins;
    return i;
}
