#include <iostream>
#include <cassert>
#include <iterator>
#include <fstream>
#include <string>
#include <algorithm>
#include <array>
#include <cmath>
#include <vector>
#include "random.hh"

typedef std::array<double, 3> point;

double squared_norm(const point & p)
{
    double s(0);
    for(double x : p)
        s += std::pow(x, 2);
    return s;
}

point operator+(const point & a, const point & b)
{
    // We should check the size of the containers...
    point c;
    for(unsigned int i = 0; i < a.size(); ++i)
        c[i] = a[i] + b[i];
    return c;
}

point operator/(const point & a, double b)
{
    point c;
    assert(b != 0 && "Division by zero.");
    for(unsigned int i = 0; i < a.size(); ++i)
        c[i] = a[i] / b;
    return c;
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

    // Start of exercise 02.1

    // 1. Random walk on a cubic lattice
    // =================================
    // The algorithm goes like this:
    // - start from the origin;
    // - at each step, move in one of the six possible directions, each one
    //   with equal probability.

    // In order to choose which way the walk will continue, we draw an
    // uniformly distributed number p on [0, 3), and
    // - if p < 1, walk along the x direction;
    // - if 1 <= p < 2, walk along the y direction;
    // - if 2 <= p < 3, walk along the z direction.
    // Then draw another number on [0, 1), and if p < 0.5 then follow the
    // "positive" direction, otherwise the "negative" one (to be defined).

    const unsigned int n_walks(10000),
                       n_steps(100);

    auto discrete_walk = [&rnd](const point & start)
    {
        // Returns the next position of the walk.
        double p;
        int axis, direction;
        p = rnd.Rannyu(0, 3);
        if(p < 1)
            axis = 0;
        if(1 <= p && p < 2)
            axis = 1;
        if(2 <= p && p < 3)
            axis = 2;
        p = rnd.Rannyu();
        if(p < 0.5)
            direction = 0;
        else
            direction = 1;
        point step(start);
        step[axis] = start[axis] + std::pow(-1, direction);
        return step;
    };
    point previous_position, position;

    std::vector<std::vector<double>> square_distance(n_walks);
    for(unsigned int i = 0; i < square_distance.size(); ++i)
        square_distance[i].resize(n_steps);

    std::vector<double> root_avg_squared_distance(n_steps);

    for(unsigned int w = 0; w < square_distance.size(); ++w)
    {
        square_distance[w][0] = 0;
        // Start from the origin.
        std::fill(previous_position.begin(), previous_position.end(), 0);
        for(unsigned int i = 1; i < root_avg_squared_distance.size(); ++i)
        {
            position = discrete_walk(previous_position);
            square_distance[w][i] = squared_norm(position);
            previous_position = position;
        }
    }

    // For each random walk we can construct an array containing, at each
    // step, the square of the distance from the starting point.
    // We can then construct the average value, taking the average over all
    // the walks, of that distance.
    double sum;
    for(unsigned int step = 0; step < root_avg_squared_distance.size(); ++step)
    {
        sum = 0; 
        for(unsigned int w = 0; w < square_distance.size(); ++w)
            sum += square_distance[w][step];
        root_avg_squared_distance[step] = std::sqrt(sum / n_walks);
    }

    std::ofstream output_file("cubic-random-walk.dat");
    std::ostream_iterator<double> out(output_file, "\n");
    std::copy(root_avg_squared_distance.begin(), root_avg_squared_distance.end(), out);

    output_file.close();
    
    // Continuous random walk
    // ======================
    // Same as above; we just need to change the walk function.

    auto uniform_walk = [&rnd](const point & start)
    {
        // Generate a point in the 2-sphere, then add it to start.
        // In order to generate the point, we use the rejection method
        // to draw a point uniformly distributed in the cube [-1, 1]^3;
        // if its norm is less than 1, we accept it.
        // This gives us a uniform distribution on the unit ball B^3:
        // normalizing the points we get a uniform distribution on S^2.
        // (We will also reject the point (0, 0, 0) in order to avoid
        // dividing by zero).
        point x;
        do
            for(double & coordinate : x)
                coordinate = rnd.Rannyu(-1, 1);
        while(squared_norm(x) > 1 || squared_norm(x) == 0);
        point step(x / std::sqrt(squared_norm(x)));

        // An alternative method, using the inversion of the cdf.
        //double theta = rnd.Rannyu(0, 2 * M_PI),
        //       phi = std::acos(1 - 2 * rnd.Rannyu());
        //point step {
        //    std::sin(phi) * std::cos(theta),
        //    std::sin(phi) * std::sin(theta),
        //    std::cos(phi)
        //};
        return start + step;
    };

    for(unsigned int w = 0; w < square_distance.size(); ++w)
    {
        square_distance[w][0] = 0;
        // Start from the origin.
        std::fill(previous_position.begin(), previous_position.end(), 0);
        for(unsigned int i = 1; i < root_avg_squared_distance.size(); ++i)
        {
            position = uniform_walk(previous_position);
            square_distance[w][i] = squared_norm(position);
            previous_position = position;
        }
    }

    // Calculation of the square of the distance from the starting point.
    for(unsigned int step = 0; step < root_avg_squared_distance.size(); ++step)
    {
        sum = 0; 
        for(unsigned int w = 0; w < square_distance.size(); ++w)
            sum += square_distance[w][step];
        root_avg_squared_distance[step] = std::sqrt(sum / n_walks);
    }

    output_file.open("continuous-random-walk.dat");
    std::copy(root_avg_squared_distance.begin(), root_avg_squared_distance.end(), out);

    rnd.SaveSeed();
    return 0;
}
