#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include "random.hh"
#include "statistics.hh"

int main(int argc, char *argv[])
{
    // Prepare the random number generator
    Random rnd;
    int seed[4];
    int p1, p2;
    std::ifstream Primes("Primes");
    if(Primes.is_open())
    {
        Primes >> p1 >> p2 ;
    }
    else
        std::cerr << "Unable to open Primes." << std::endl;
    Primes.close();

    std::ifstream input("seed.in");
    std::string property;
    if(input.is_open())
    {
        while(!input.eof())
        {
            input >> property;
            if(property == "RANDOMSEED")
            {
                input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
                rnd.SetRandom(seed,p1,p2);
            }
        }
        input.close();
    }
    else
        std::cerr << "Unable to open seed.in." << std::endl;

    // Start of exercise 02.1

    /*
       1. Basic Monte Carlo integration
       ================================
       - draw an uniformely distributed number on the domain of the function
       to be integrated, g;
       - compute g(x) of such number x;
       - the average value of all such g(x)s is an estimate of the integral.

       2. Monte Carlo integration with importance sampling
       ===================================================
       The integration algorithm is as the one before, except that the number
       is drawn from a distribution of choice, f, and calculating the average
       of g(x) / f(x) instead of g(x) alone.
       A good sampling function for the rejection method would be the
       Taylor series of the integrand function g, truncated to o(x^3),
       with minor corrections: we start from
       f(x) = a * (1 - (x * pi / 2)^2 / 2 + b),
       with a and b to be determined such that f is as close as possible
       to the integrand. An optimal choice is such that f(1) = g(1) = 0;
       this and the normalisation to 1 gives
       a = 12 / pi^2,
       b = pi^2 / 8 - 1.
       The sampling pdf is then f(x) = 3 / 2 * (1 - x^2).
       Its cdf is F(x) = 3 / 2 * (x - x^3 / 3); inverting this function on the
       interval [0, 1] gives us the function
       p(x) = 2 * cos(5 * pi / 3 - 1 / 3 * arccos(x)),
       which we will use to generate the samples.
    */

    const int block_size(1e3),
              n_draws(1e5);
    
    std::vector<double> standard_samples,
                        standard_integral_avg,
                        standard_integral_std,
                        importance_samples,
                        importance_integral_avg,
                        importance_integral_std;
    
    standard_samples.reserve(n_draws);
    importance_samples.reserve(n_draws);

    auto integrand_function = [](double x) {return M_PI / 2 * std::cos(M_PI * x / 2);};

    auto pdf = [](double x)
    {
        return 2. * (1 - x);
    };
    auto draw = [&rnd]()
    {
        // Returns the inverse of the cumulative function of the pdf above.
        return 1. + std::sqrt(1 - rnd.Rannyu());
    };

    double x;
    for(unsigned int i = 0; i < n_draws; ++i)
    {
        standard_samples.push_back(integrand_function(rnd.Rannyu()));
        x = draw();
        importance_samples.push_back(integrand_function(x) / pdf(x));
    }

    block_statistics(
            std::begin(standard_samples),
            std::end(standard_samples),
            std::back_inserter(standard_integral_avg),
            std::back_inserter(standard_integral_std),
            block_size
            );

    block_statistics(
            std::begin(importance_samples),
            std::end(importance_samples),
            std::back_inserter(importance_integral_avg),
            std::back_inserter(importance_integral_std),
            block_size
            );

    std::ofstream output_file("integral_standard_sampling.dat");
    const unsigned int col_width = 16;

    // Header row
    output_file << std::setw(col_width) << "integral_avg" << std::setw(col_width) << "integral_std" << std::endl;

    output_file.precision(4);
    output_file << std::scientific;

    for(unsigned int i = 0; i < standard_integral_avg.size(); ++i)
        output_file << std::setw(col_width) << standard_integral_avg[i] << std::setw(col_width) << standard_integral_std[i] << std::endl;

    output_file.close();

    output_file.open("integral_importance_sampling.dat");

    // Header row
    output_file << std::setw(col_width) << "integral_avg" << std::setw(col_width) << "integral_std" << std::endl;

    output_file.precision(4);
    output_file << std::scientific;

    for(unsigned int i = 0; i < importance_integral_avg.size(); ++i)
        output_file << std::setw(col_width) << importance_integral_avg[i] << std::setw(col_width) << importance_integral_std[i] << std::endl;

    output_file.close();

    rnd.SaveSeed();
    return 0;
}
