#ifndef RANDOM_HH
#define RANDOM_HH

#include <string>

namespace lsn
{

class random
{
	public:
		random(const std::string & primer_file, const std::string & seed_file);
		~random();

		void save_seed() const;

		// Common distributions
		double uniform();
		double uniform(double min, double max);
		double normal(double mean, double stdev);

	private:
		int m1, m2, m3, m4,
			l1, l2, l3, l4,
			n1, n2, n3, n4;
};

} // namespace lsn
#endif
