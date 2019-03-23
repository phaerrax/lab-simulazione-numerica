#include <iostream>
#include <fstream>
#include <string>
#include <numeric>
#include <cmath>
#include <vector>
#include "random.hh"

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

    std::vector<double> r(total_n_throws);
    // The test
    // --------
    // We draw a certain number of (supposedly) uniformly distributed numbers.
    // For each draw, we save it in a vector, we calculate its squared
    // distance from the expected mean value (i.e. 1/2); from these data we
    // can compute the average value and the standard deviation of the set.
    for(auto p = r.begin(); p != r.end(); ++p)
        *p = rnd.Rannyu();

    std::vector<int> x(n_blocks); // Block indices.
    // Fill x with sequentially increasing values, starting from 0 and
    // adding 1 each time.
    std::iota(x.begin(), x.end(), 0);

    std::vector<double> average(n_blocks), // Average value WITHIN each block.
                        variance(n_blocks), // Variance WITHIN each block.
                        average_avg(n_blocks), // Average of all the blocks.
                        variance_avg(n_blocks),
                        average_std(n_blocks), // St dev of all the blocks.
                        variance_std(n_blocks),
                        squared_average(n_blocks),
                        squared_variance(n_blocks);

    double partial_sum_average,
           partial_sum_variance;
    for(unsigned int i = 0; i < x.size(); ++i) // Loop for each block.
    {
        // Compute the average of the drawn numbers in the i-th block.
        partial_sum_average  = 0;
        partial_sum_variance = 0;
        for(unsigned int j = 0; j < n_throws; ++j)
        {
            partial_sum_average  += r[n_throws * i + j];
            partial_sum_variance += std::pow(r[n_throws * i + j] - 0.5, 2);
        }
        average[i]          = partial_sum_average / n_throws; 
        variance[i]         = partial_sum_variance / n_throws;
        // For each block the std deviation of the average value from 1/2
        // can be computed using the formula given at the end of the notebook;
        // ditto for the variance, whose expected value is 1/12.
        // For this purpose we will need the squares of the values above.
        squared_average[i]  = std::pow(average[i], 2);
        squared_variance[i] = std::pow(variance[i], 2);
    }

    double sum_average, // Average of the avg values of the first i blocks.
           sum_variance, // Same but for the variance.
           sum_squared_average,
           sum_squared_variance;
    for(unsigned int i = 0; i < x.size(); ++i)
    {
        sum_average = 0;
        sum_variance = 0;
        sum_squared_average = 0;
        sum_squared_variance = 0;
        for(unsigned int j = 0; j <= i; ++j)
        {
            sum_average          += average[j];
            sum_variance         += variance[j];
            sum_squared_average  += std::pow(average[j], 2);
            sum_squared_variance += std::pow(variance[j], 2);
        }
        average_avg[i]  = sum_average / (i + 1); // (i runs from 0 to n_blocks - 1.)
        variance_avg[i] = sum_variance / (i + 1);
        if(i > 0)
        {
            average_std[i]  = std::sqrt((sum_squared_average / (i + 1) - std::pow(average_avg[i], 2)) / i);
            variance_std[i] = std::sqrt((sum_squared_variance / (i + 1)- std::pow(variance_avg[i], 2)) / i);
        }
        else
        {
            average_std[i]  = 0;
            variance_std[i] = 0;
        }
    }

    for(unsigned int i = 0; i < x.size(); ++i)
        x[i] *= n_throws;

    // Output
    std::ofstream average_output("uniformd-average.dat"),
                  variance_output("uniformd-variance.dat");
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
    std::ofstream output("uniformd-chisquared-test.dat");
    for(unsigned int i = 0; i < chisquared_values.size(); ++i)
        output << i << " " << chisquared_values[i] << std::endl;
    output.close();

    rnd.SaveSeed();
    return 0;
}
