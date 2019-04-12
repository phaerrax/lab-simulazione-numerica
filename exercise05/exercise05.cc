#include <cmath>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <iterator>
#include <numeric>
#include <functional>
#include <array>
#include <vector>
#include "metropolis.hh"
#include "metropolis_uniform.hh"
#include "random.hh"

std::vector<std::vector<double>> block_statistics(const std::vector<double> &, unsigned int);

int main()
{
    // Random number generator initialization
    Random rng;
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
                rng.SetRandom(seed,p1,p2);
            }
        }
        input.close();
    }
    else
        std::cerr << "Unable to open seed.in." << std::endl;

    const unsigned int dim = 3;
    double L(1.2);
    metropolis_uniform<dim> metro(-L, L);
    
    std::array<double, dim> start = {1, 0, 0};
    metro.set_starting_point(start);

    auto radius = [](const std::array<double, dim> & x)
    {
        return std::sqrt(
                std::inner_product(x.begin(), x.end(), x.begin(), 0.)
                );
    };

    std::function<double (std::array<double, dim>)> f = [radius](std::array<double, dim> x)
    {
        return std::exp(-2. * radius(x)) / M_PI;
    };
    // The function to be sampled.

    unsigned int n_steps(1e5);
    std::vector<std::array<double, dim>> sequence;
    std::array<double, dim> new_point;
    std::vector<double> r;
    for(unsigned int n = 0; n < n_steps; ++n)
    {
        new_point = metro.step(f, rng);
        sequence.push_back(new_point);
        r.push_back(radius(new_point));
    }

    std::cerr << "Acceptance rate after " << n_steps << " steps: " << metro.get_acceptance_rate() << "." << std::endl;

    // Output the sequence of sampled points to a file.
    std::ofstream output_file("sampled_points_uniform.dat");

    const unsigned int col_width = 16;
    output_file.precision(4);
    output_file << std::scientific;

    for(auto row : sequence)
    {
        for(auto x = row.begin(); x != row.end(); ++x)
            output_file << std::setw(col_width) << *x;
        output_file << "\n";
    }
    output_file << std::endl;

    output_file.close();

    // Output the progressive values of the average radius and its standard
    // deviation, obtained with a blocking technique, to a file.
    output_file.open("radius_avg_uniform.dat");

    std::vector<std::vector<double>> r_blocks(block_statistics(r, 100));

    for(auto row : r_blocks)
        output_file << std::setw(col_width) << row[0]
                    << std::setw(col_width) << row[1]
                    << "\n";
    output_file << std::endl;

    output_file.close();

    return 0;
}

std::vector<std::vector<double>> block_statistics(const std::vector<double> & x, unsigned int n_blocks)
{
    unsigned int block_size = static_cast<unsigned int>(std::round(
            static_cast<double>(x.size()) / n_blocks
            ));

    // Just to make sure:
    assert(block_size * n_blocks == x.size());

    double sum(0), sum_sq(0), block_average;
    std::vector<double> row(2);
    std::vector<std::vector<double>> result;

    // - sum the values in each block;
    // - compute the average of that block;
    // - from the list of averages compute the standard dev of the mean.
    for(unsigned int i = 0; i < n_blocks; ++i)
    {
        block_average = 0;
        for(unsigned int j = 0; j < block_size; ++j)
            block_average += x[i * block_size + j];
        block_average /= block_size;
        sum           += block_average;
        sum_sq        += std::pow(block_average, 2);
        row[0]         = sum / (i + 1);
        if(i > 0)
            row[1] = std::sqrt((sum_sq / (i + 1) - std::pow(row[0], 2)) / i);
        else
            row[1] = 0;
        result.push_back(row);
    }
    return result;
}
