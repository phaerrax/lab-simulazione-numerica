#include <iostream>
#include <fstream>
#include <string>
#include <numeric>
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

    // Start of exercise 01.1

    const int n_blocks(100);
    const int n_throws(1000); // (Inside each block.)
    int total_n_throws = n_throws * n_blocks;

    std::vector<double> r(total_n_throws),
                        r_var;
    // The test
    // --------
    // We draw a certain number of (supposedly) uniformly distributed numbers.
    // For each draw, we save it in a vector, we calculate its squared
    // distance from the expected mean value (i.e. 1/2); from these data we
    // can compute the average value and the standard deviation of the set.
    for(auto & p : r)
    {
        p = rnd.Rannyu();
        r_var.push_back(std::pow(p - 0.5, 2));
    }

    std::vector<int> x(n_blocks); // Block indices.
    // Fill x with sequentially increasing values, starting from 0 and
    // adding 1 each time.
    std::iota(x.begin(), x.end(), 0);

    std::vector<double> average_avg, // Average of all the blocks.
                        variance_avg,
                        average_std, // St dev of all the blocks.
                        variance_std;

    block_statistics(
            std::begin(r),
            std::end(r),
            std::back_inserter(average_avg),
            std::back_inserter(average_std),
            n_throws
            );

    block_statistics(
            std::begin(r_var),
            std::end(r_var),
            std::back_inserter(variance_avg),
            std::back_inserter(variance_std),
            n_throws
            );

    for(unsigned int i = 0; i < x.size(); ++i)
        x[i] *= n_throws;

    // Output
    std::ofstream average_output("uniformd_average.dat"),
                  variance_output("uniformd_variance.dat");
    for(unsigned int i = 0; i < x.size(); ++i)
    {
        average_output << x[i] << " " << average_avg[i] << " " << average_std[i] << std::endl;
        variance_output << x[i] << " " << variance_avg[i] << " " << variance_std[i] << std::endl;
    }
    average_output.close();
    variance_output.close();

    // Start of exercise 01.1 - part 3

    // Procedure:
    // - if the distribution is truly uniform, each subinterval should
    //   contain N/M numbers on average;
    // - we compare in each bin (i = 0, ..., M-1) how many numbers
    //   drawn previously fall in such bin (this is the observed
    //   value) with the expected value N/M;
    // - we calculate then the chi-squared for the experiment.

    const int n_bins(100),
              n_tries = n_blocks,
              n_draws = n_throws;
    // These assingments are vital if we plan to test the numbers in the
    // vector r that we previously defined, which has length exactly
    // n_blocks * n_throws.

    double sum, draw;

    std::vector<double> chisquared_values(n_tries);
    std::vector<int> bin_count(n_bins);

    double expected_value = static_cast<double>(n_draws)/n_bins;

    auto p = r.begin();

    for(int i = 0; i < n_tries; ++i)
    {
        sum = 0;
        // Reset bin_count.
        std::fill(bin_count.begin(), bin_count.end(), 0);
        for(int j = 0; j < n_draws; ++j)
        {
            // Get the next uniformly drawn number from the list.
            draw = *p;
            ++p;
            // If we multiply by n_bins and take the integer part, we should
            // get the number of the bin in which the number is located.
            bin_count[static_cast<int>(std::floor(draw * n_bins))]++;
        }
        // Compute chi-squared from the bin count.
        for(unsigned int j = 0; j < bin_count.size(); ++j)
            sum += std::pow(static_cast<double>(bin_count[j]) - expected_value, 2);
        chisquared_values[i] = sum / expected_value;
    }

    // Output
    std::ofstream output("uniformd_chisquared_test.dat");
    for(unsigned int i = 0; i < chisquared_values.size(); ++i)
        output << i << " " << chisquared_values[i] << std::endl;
    output.close();

    rnd.SaveSeed();
    return 0;
}
