/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <iostream>
#include <fstream>
#include <cmath>
#include "random.hh"

namespace lsn
{

random::random(const std::string & primer_file, const std::string & seed_file)
{
	int seed[4];
	int p1, p2;
	std::ifstream Primes(primer_file);
	if(Primes.is_open())
	{
		Primes >> p1 >> p2 ;
	}
	else
		std::cerr << "Unable to open Primes." << std::endl;
	Primes.close();

	std::ifstream input(seed_file);
	std::string property;
	if(input.is_open())
	{
		while(!input.eof())
		{
			input >> property;
			if(property == "RANDOMSEED")
			{
				input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
				m1 = 502;
				m2 = 1521;
				m3 = 4071;
				m4 = 2107;
				l1 = seed[0] % 4096;
				l2 = seed[1] % 4096;
				l3 = seed[2] % 4096;
				l4 = seed[3] % 4096;
				l4 = 2 * (l4 / 2) + 1;
				n1 = 0;
				n2 = 0;
				n3 = p1;
				n4 = p2;
			}
		}
		input.close();
	}
	else
		std::cerr << "Unable to open seed.in." << std::endl;
}

random::~random()
{
	save_seed();
}

void random::save_seed() const
{
	std::ofstream write_seed;
	const std::string output_file("seed.out");
	write_seed.open(output_file);
	if(write_seed.is_open())
	{
		write_seed << l1 << " " << l2 << " " << l3 << " " << l4 << std::endl;
	}
	else
		std::cerr << "Unable to open " << output_file << "." << std::endl;
	write_seed.close();
	return;
}

double random::normal(double mean, double stdev)
{
	double s(uniform()),
	       t(uniform());
	double x = std::sqrt(-2. * std::log(1. - s)) * std::cos(2. * M_PI * t);
	return mean + x * stdev;
}

double random::uniform(double min, double max)
{
	return min + (max - min) * uniform();
}

double random::uniform()
{
	const double twom12(0.000244140625);
	int i1, i2, i3, i4;
	double r;

	i1 = l1*m4 + l2*m3 + l3*m2 + l4*m1 + n1;
	i2 = l2*m4 + l3*m3 + l4*m2 + n2;
	i3 = l3*m4 + l4*m3 + n3;
	i4 = l4*m4 + n4;
	l4 = i4 % 4096;
	i3 = i3 + i4 / 4096;
	l3 = i3 % 4096;
	i2 = i2 + i3 / 4096;
	l2 = i2 % 4096;
	l1 = (i1 + i2 / 4096) % 4096;
	r = twom12 * (l1 + twom12 * (l2 + twom12 * (l3 + twom12 * l4)));

	return r;
}

} // namespace lsn
