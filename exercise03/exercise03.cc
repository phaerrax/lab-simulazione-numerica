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

    // Start of exercise 03.1
	
	// 1. Compute the price of a call or put option price by a direct sampling
	// of the price at the expiry.
	// =======================================================================
	//
	// - sample N times a geometric Brownian motion s(T) at the expiry T;
	// - from the sample, estimate the option price as the profit
	//   max(S(T)-K,0), K being the strike price;
	// - repeat a certain number of times, each time obtaining a new estimate;
	// - calculate the progression of the average value and its standard
	//   deviation as the number of blocks increases.
	
	const double expiry(1),
				 starting_asset_price(100),
				 strike_price(100),
				 risk_free_interest_rate(0.1),
				 volatility(0.25);

	const unsigned int n_data(1e6),
					   block_size(1e4);
	
    std::vector<double> call_profit,
                        put_profit,
                        avg_call_option_price,
	                    std_call_option_price,
	                    avg_put_option_price,
	                    std_put_option_price;

    call_profit.reserve(n_data);
    put_profit.reserve(n_data);

	double expiry_price;
	for(unsigned int i = 0; i < n_data; ++i)
	{
        // Sample the final asset price with a geometric Browian motion
        // simulation.
        expiry_price = starting_asset_price * std::exp(expiry * (risk_free_interest_rate - 0.5 * std::pow(volatility, 2)) + volatility * rnd.Gauss(0, expiry));
        call_profit.push_back(std::exp(-risk_free_interest_rate * expiry) * std::max(expiry_price - strike_price, 0.));
        put_profit.push_back(std::exp(-risk_free_interest_rate * expiry) * std::max(strike_price - expiry_price, 0.));
    }

    block_statistics(
            std::begin(call_profit),
            std::end(call_profit),
            std::back_inserter(avg_call_option_price),
            std::back_inserter(std_call_option_price),
            block_size
            );

    block_statistics(
            std::begin(put_profit),
            std::end(put_profit),
            std::back_inserter(avg_put_option_price),
            std::back_inserter(std_put_option_price),
            block_size
            );

	std::ofstream call_output_file("call_direct_price.dat"),
	              put_output_file("put_direct_price.dat");
	for(unsigned int i = 0; i < avg_call_option_price.size(); ++i)
	{
		call_output_file << avg_call_option_price[i] << " " << std_call_option_price[i] << "\n";
		put_output_file  << avg_put_option_price[i]  << " " << std_put_option_price[i]  << "\n";
	}
	call_output_file << std::endl;
	put_output_file  << std::endl;

	call_output_file.close();
	put_output_file.close();

	// 2. Compute the price of a call or put option price by sampling
	//    progressively a geometric Brownian motion.
	// =======================================================================
	//
	// - sample N sequences of a geometric Brownian motion s(t_i) from t_0 = 0
	//   to t_n = N;
	// - the rest goes on as above.

	unsigned int n_steps(100); // Steps to sample the random path of the price.
	double time_step(expiry / n_steps);

    call_profit.clear();
    put_profit.clear();
    avg_call_option_price.clear();
    std_call_option_price.clear();
    avg_put_option_price.clear();
    std_put_option_price.clear();

	for(unsigned int i = 0; i < n_data; ++i)
	{
        expiry_price = starting_asset_price;
        // Sample the asset price at the expiry by sampling the process
        // along its random path, with small increments.
        for(unsigned int step = 0; step < n_steps; ++step)
            expiry_price *= std::exp((risk_free_interest_rate - 0.5 * std::pow(volatility, 2)) * time_step + volatility * rnd.Gauss(0, 1) * std::sqrt(time_step));
        call_profit.push_back(std::exp(-risk_free_interest_rate * expiry) * std::max(expiry_price - strike_price, 0.));
        put_profit.push_back(std::exp(-risk_free_interest_rate * expiry) * std::max(strike_price - expiry_price, 0.));
    }

    block_statistics(
            std::begin(call_profit),
            std::end(call_profit),
            std::back_inserter(avg_call_option_price),
            std::back_inserter(std_call_option_price),
            block_size
            );

    block_statistics(
            std::begin(put_profit),
            std::end(put_profit),
            std::back_inserter(avg_put_option_price),
            std::back_inserter(std_put_option_price),
            block_size
            );

	call_output_file.open("call_progressive_price.dat"),
	put_output_file.open("put_progressive_price.dat");
	for(unsigned int i = 0; i < avg_call_option_price.size(); ++i)
	{
		call_output_file << avg_call_option_price[i] << " " << std_call_option_price[i] << "\n";
		put_output_file  << avg_put_option_price[i]  << " " << std_put_option_price[i]  << "\n";
	}
	call_output_file << std::endl;
	put_output_file  << std::endl;

	call_output_file.close();
	put_output_file.close();

    rnd.SaveSeed();
    return 0;
}
