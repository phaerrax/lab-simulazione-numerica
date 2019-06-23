#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include "random.hh"
#include <iomanip>

std::string filename(std::string dist_name)
{
    // Standard filename scheme of the output files.
    return "central_limit_theorem_" + dist_name + "d.dat";
}

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

    // Start of exercise 01.2

    // Procedure:
    // - draw M numbers following the distribution X;
    // - compute their average;
    // - repeat T times;
    // - plot the distribution of the set of T averages;
    // - maybe also vary M and see what changes.

    // Uniform distribution
    const int n_tries = 10000; // "number of realizations of S_N"
    const std::vector<unsigned int> n_draws {1, 2, 10, 100};
    const std::vector<std::string> dist_name {"uniform", "exponential", "cauchylorentz"};
    std::vector<std::vector<std::vector<double>>> averages(dist_name.size());
    // (0) uniform distribution
    // (1) exponential distribution
    // (2) Cauchy-Lorentz distribution

    // Each averages[i] will be a matrix whose columns are indexed by n_draws; each element
    // in the row N is an average of N elements following the associated distribution.
    for(unsigned int i = 0; i < dist_name.size(); ++i)
    {
        // Rows initializazion.
        averages[i].resize(n_draws.size());
        // Columns inizialization.
        for(unsigned int j = 0; j < averages[i].size(); ++j)
            averages[i][j].resize(n_tries, 0);
    }

    unsigned int N;
    for(unsigned int i = 0; i < n_draws.size(); ++i)
    {
        N = n_draws[i];
        for(unsigned int k = 0; k < n_tries; ++k)
        {
            for(unsigned int j = 0; j < N; ++j)
            {
                averages[0][i][k] += rnd.Rannyu();
                averages[1][i][k] += rnd.exponential(1);
                averages[2][i][k] += rnd.cauchylorentz(0, 1);
            }
            for(unsigned int type = 0; type < dist_name.size(); ++type)
                averages[type][i][k] /= N;
        }
    }

    // Output
    std::vector<std::ofstream> output_files(dist_name.size());
    for(unsigned int i = 0; i < dist_name.size(); ++i)
        output_files[i].open(filename(dist_name[i]));

    const unsigned int col_width = 12;

    // Header row
    for(unsigned int i = 0; i < dist_name.size(); ++i)
    {
        for(unsigned int n = 0; n < n_draws.size(); ++n)
            output_files[i] << std::setw(col_width) << n_draws[n];
        output_files[i] << std::endl;
    }


    for(unsigned int dist_type = 0; dist_type < dist_name.size(); ++dist_type)
    {
        output_files[dist_type].precision(4);
        output_files[dist_type] << std::scientific;
        for(unsigned int k = 0; k < n_tries; ++k)
        {
            for(unsigned int i = 0; i < n_draws.size(); ++i)
                output_files[dist_type] << std::setw(col_width) << averages[dist_type][i][k];
            output_files[dist_type] << std::endl;
        }
    }

    for(unsigned int i = 0; i < dist_name.size(); ++i)
        output_files[i].close();

    rnd.SaveSeed();
    return 0;
}
