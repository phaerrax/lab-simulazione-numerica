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
#include "statistics.hh"

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

    for(const auto & row : sequence1s)
    {
        for(auto x = row.begin(); x != row.end(); ++x)
            output_file << std::setw(col_width) << *x;
        output_file << "\n";
    }
    output_file << std::endl;
    output_file.close();

    output_file.open("2p_sampled_points_uniform.dat");
    for(const auto & row : sequence2p)
    {
        for(auto x = row.begin(); x != row.end(); ++x)
            output_file << std::setw(col_width) << *x;
        output_file << "\n";
    }
    output_file << std::endl;
    output_file.close();

    // Output the progressive values of the average radius and its standard
    // deviation, obtained with a blocking technique, to a file.
    std::vector<double> r1s_avg,
                        r1s_std,
                        r2p_avg,
                        r2p_std;

    output_file.open("1s_radius_avg_uniform.dat");
    block_statistics(
            std::begin(r1s),
            std::begin(r1s),
            std::back_inserter(r1s_avg),
            std::back_inserter(r1s_std),
            r1s.size() / 100
            );
    for(unsigned int i = 0; i < r1s.size(); ++i)
        output_file << std::setw(col_width) << r1s_avg[i]
                    << std::setw(col_width) << r1s_std[i]
                    << "\n";
    output_file << std::endl;
    output_file.close();

    output_file.open("2p_radius_avg_uniform.dat");
    block_statistics(
            std::begin(r2p),
            std::begin(r2p),
            std::back_inserter(r2p_avg),
            std::back_inserter(r2p_std),
            r2p.size() / 100
            );
    for(unsigned int i = 0; i < r2p.size(); ++i)
        output_file << std::setw(col_width) << r2p_avg[i]
                    << std::setw(col_width) << r2p_std[i]
                    << "\n";
    output_file << std::endl;
    output_file.close();

    return 0;
}
