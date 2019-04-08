#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <vector>
#include "random.hh"
#include "molecular_dynamics_sim.hh"

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

    Input();             //Inizialization
    int nconf = 1;
    for(int istep=1; istep <= nstep; ++istep)
    {
        Move();           //Move particles with Verlet algorithm
        if(istep%iprint == 0)
            std::cout << "Number of time-steps: " << istep << std::endl;
        if(istep%10 == 0)
        {
            Measure();     //Properties measurement
            //        ConfXYZ(nconf);//Write actual configuration in XYZ format //Commented to avoid "filesystem full"! 
            nconf += 1;
        }
    }
    ConfFinal();         //Write final configuration to restart

    return 0;
}
