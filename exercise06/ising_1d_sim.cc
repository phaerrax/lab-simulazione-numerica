#include "ising_1d_sim.hh"
#include "random.hh"
#include <numeric>
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
	double acceptance_threshold,
		   energy_diff;
	for(unsigned int i = 0; i < spins.size(); ++i)
	{
		energy_diff = 2 * spins[quotient(i)] * (interaction_energy * (spins[quotient(i + 1)] + spins[quotient(i - 1)]) + ext_magnetic_field);
		acceptance_threshold = std::exp(-inverse_temperature * energy_diff);
		if(acceptance_threshold < 1)
		{
			if(rng.Rannyu() < acceptance_threshold)
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
		// The probability that the spin will be -1 is 1 / (1 + e^(-2 * x)).
		if(rng.Rannyu() * (1. + std::exp(-2. * x)) < 1.)
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
