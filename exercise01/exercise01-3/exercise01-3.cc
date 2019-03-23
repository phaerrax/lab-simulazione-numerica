#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <cmath>
#include <complex>
#include <vector>
#include <iterator>
#include <algorithm>
#include <utility>
#include <array>
#include "random.hh"

typedef std::complex<double> point;

std::pair<point, point> throw_stick(Random &, double);
point uniformd_circle(Random &, point, double);

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

    // Start of exercise 01.3

    const double L = 0.5; // The length of the stick.

    const unsigned int n_tries(100), n_throws(10000);
    std::vector<double> pi_estimate(n_tries),     // Average within a block.
                        pi_estimate_avg(n_tries), // Overall average.
                        pi_estimate_std(n_tries); // Overall std.
    unsigned int hits;
    double sum_averages(0), sum_squared_averages(0);
    std::pair<point, point> p;

    // Make M throws and for each step compute the average and the
    // uncertainty.
    // Repeat N times: as more and more results are added, the average
    // should improve since the (overall) uncertainty tends to zero.
    for(unsigned int i = 0; i < n_tries; ++i)
    {
        hits = 0;
        for(unsigned int j = 0; j < n_throws; ++j)
        {
            p = throw_stick(rnd, L);
            // We test if the stick crosses a line.
            // Since the lines are at x = 0, 1, 2, ..., 10 this happens if the
            // integer parts of the x-coordinates of the points are different.
            if(std::floor(std::real(p.first)) != std::floor(std::real(p.second)))
                ++hits;
        }
        pi_estimate[i] = 2 * L * n_throws / hits; // (d = 1.)
        // Maybe we can write here a check in case hits turns out to be zero...

        // At each step i, we compute the sum of the average values (i.e. the
        // estimates of pi) and of their squares. We will need these to
        // calculate the standard deviation at each step.
        sum_averages += pi_estimate[i];
        sum_squared_averages += std::pow(pi_estimate[i], 2);

        pi_estimate_avg[i] = sum_averages / (i + 1);
        // (i runs from 0 to n_blocks - 1.)
        if(i > 0)
            pi_estimate_std[i] = std::sqrt((sum_squared_averages / (i + 1) - std::pow(pi_estimate_avg[i], 2)) / i);
        else
            pi_estimate_std[i] = 0;
    }

    std::ofstream output_file("pi_estimate.dat");
    const unsigned int col_width = 12;

    // Header row
    output_file << std::setw(col_width) << "pi_estimate" << std::setw(col_width) << "st_dev" << std::endl;

    output_file.precision(4);
    output_file << std::scientific;

    for(unsigned int i = 0; i < n_tries; ++i)
        output_file << std::setw(col_width) << pi_estimate_avg[i] << std::setw(col_width) << pi_estimate_std[i] << std::endl;

    output_file.close();

    rnd.SaveSeed();
    return 0;
}

point uniformd_circle(Random & rng, point centre, double radius)
{
    double t1, t2;
    do
    {
        t1 = rng.Rannyu(-1, 1);
        t2 = rng.Rannyu(-1, 1);
    }
    while(std::hypot(t1, t2) > 1);
    // If the point lies outside the circle, draw another one.
    // If it is inside, we are done and can return it.
    point t(std::real(centre) + t1, std::imag(centre) + t2);
    return t;
}

std::pair<point, point> throw_stick(Random & rng, double L)
{
    // (L: length of the stick)
    //
    // The floor (x,y) is ruled with a vertical line at each integer x.
    // We will draw the coordinates from a uniform distribution, so it
    // must be bounded: for this reason we consider as "the floor" the
    // square [0,10]x[0,10].
    //
    // ## How to simulate the throw of a stick, without using pi ##
    // We draw two random points on the plane and normalize the length
    // such that it is the assinged length of the stick.
    // The first one is chosen at random in our previously defined "floor";
    // the second one is drawn from a uniform distribution on a circle
    // centered on the first point and of radius 1.
    // To sample this we use a simple rejection algorithm.

    // Draw the coordinates of the first point.
    point a(rng.Rannyu(0, 10), rng.Rannyu(0, 10));
    // Draw the second from a circle around the first.
    point b;
    b = uniformd_circle(rng, a, 1);
    // Now, normalize b such that it is at a distance L from a.
    // The result is the point (1 - t) * a + t * b where t is defined
    // as L^2 / |a - b|^2.
    b = a + L / std::abs(b - a) * (b - a);
    // The stick is now a segment [a,b] of length L.

    std::pair<point, point> p(a,b);
    return p;
}
