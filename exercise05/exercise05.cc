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
    metropolis_uniform<dim> metro1s(-L, L);
    metropolis_uniform<dim> metro2p(-L, L);
    
    std::array<double, dim> start = {1, 0, 0};
    metro1s.set_starting_point(start);
    metro2p.set_starting_point(start);

    auto radius = [](const std::array<double, dim> & x)
    {
        return std::sqrt(
                std::inner_product(x.begin(), x.end(), x.begin(), 0.)
                );
    };

    std::function<double (std::array<double, dim>)> f1s = [radius](std::array<double, dim> x)
    {
        return std::exp(-2. * radius(x)) / M_PI;
    };
    std::function<double (std::array<double, dim>)> f2p = [radius](std::array<double, dim> x)
    {
        return std::pow(x[2], 2) * std::exp(-radius(x)) / (32 * M_PI);
    };
    // The functions to be sampled.

    unsigned int n_steps(1e6);
    std::vector<std::array<double, dim>> sequence1s, sequence2p;
    std::array<double, dim> new_point;
    std::vector<double> r1s, r2p;
    for(unsigned int n = 0; n < n_steps; ++n)
    {
        new_point = metro1s.step(f1s, rng);
        sequence1s.push_back(new_point);
        r1s.push_back(radius(new_point));

        new_point = metro2p.step(f2p, rng);
        sequence2p.push_back(new_point);
        r2p.push_back(radius(new_point));
    }

    std::cerr << "[1s] Acceptance rate after " << n_steps << " steps: " << metro1s.get_acceptance_rate() << "." << std::endl;
    std::cerr << "[2p] Acceptance rate after " << n_steps << " steps: " << metro2p.get_acceptance_rate() << "." << std::endl;

    // Output the sequence of sampled points to a file.
    std::ofstream output_file("1s_sampled_points_uniform.dat");

    const unsigned int col_width = 16;
    output_file.precision(4);
    output_file << std::scientific;

    for(auto row : sequence1s)
    {
        for(auto x = row.begin(); x != row.end(); ++x)
            output_file << std::setw(col_width) << *x;
        output_file << "\n";
    }
    output_file << std::endl;
    output_file.close();

    output_file.open("2p_sampled_points_uniform.dat");
    for(auto row : sequence2p)
    {
        for(auto x = row.begin(); x != row.end(); ++x)
            output_file << std::setw(col_width) << *x;
        output_file << "\n";
    }
    output_file << std::endl;
    output_file.close();

    // Output the progressive values of the average radius and its standard
    // deviation, obtained with a blocking technique, to a file.
    output_file.open("1s_radius_avg_uniform.dat");
    std::vector<std::vector<double>> r1s_blocks(block_statistics(r1s, 100));
    for(auto row : r1s_blocks)
        output_file << std::setw(col_width) << row[0]
                    << std::setw(col_width) << row[1]
                    << "\n";
    output_file << std::endl;
    output_file.close();

    output_file.open("2p_radius_avg_uniform.dat");
    std::vector<std::vector<double>> r2p_blocks(block_statistics(r2p, 100));
    for(auto row : r2p_blocks)
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
