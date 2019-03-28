#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
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

	const unsigned int n_blocks(100),
					   block_size(10000);

	double call_estimates_sum(0),
		   call_squared_estimates_sum(0),
	       put_estimates_sum(0),
		   put_squared_estimates_sum(0);
	
	std::vector<double> avg_call_option_price(n_blocks),
	                    avg_put_option_price(n_blocks),
	                    std_call_option_price(n_blocks),
	                    std_put_option_price(n_blocks);

	double call_profit, put_profit, expiry_price;
	for(unsigned int i = 0; i < n_blocks; ++i)
	{
		call_profit = 0;
		put_profit  = 0;
		for(unsigned int j = 0; j < block_size; ++j)
		{
			// Sample the final asset price with a geometric Browian motion
			// simulation.
			expiry_price = starting_asset_price * std::exp(expiry * (risk_free_interest_rate - 0.5 * std::pow(volatility, 2)) + volatility * rnd.Gauss(0, expiry));
			call_profit  += std::exp(-risk_free_interest_rate * expiry) * std::max(expiry_price - strike_price, 0.);
			put_profit   += std::exp(-risk_free_interest_rate * expiry) * std::max(strike_price - expiry_price, 0.);
		}
		call_profit /= block_size;
		call_estimates_sum += call_profit;
		call_squared_estimates_sum += std::pow(call_profit, 2);
		avg_call_option_price[i] = call_estimates_sum / (i + 1);
		put_profit /= block_size;
		put_estimates_sum += put_profit;
		put_squared_estimates_sum += std::pow(put_profit, 2);
		avg_put_option_price[i] = put_estimates_sum / (i + 1);
		if(i > 0)
		{
			std_call_option_price[i] = std::sqrt((call_squared_estimates_sum / (i + 1) - std::pow(avg_call_option_price[i], 2)) / i);
			std_put_option_price[i]  = std::sqrt((put_squared_estimates_sum / (i + 1) - std::pow(avg_put_option_price[i], 2)) / i);
		}
		else
		{
			std_call_option_price[i] = 0;
			std_put_option_price[i]  = 0;
		}
	}

	std::ofstream call_output_file("call-direct-price.dat"),
	              put_output_file("put-direct-price.dat");
	for(unsigned int i = 0; i < n_blocks; ++i)
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

	call_estimates_sum = 0;
	call_squared_estimates_sum = 0;
	put_estimates_sum = 0;
	put_squared_estimates_sum = 0;
	
	unsigned int n_steps(100); // Steps to sample the random path of the price.
	double time_step(expiry / n_steps);

	for(unsigned int i = 0; i < n_blocks; ++i)
	{
		call_profit = 0;
		put_profit  = 0;
		for(unsigned int j = 0; j < block_size; ++j)
		{
			expiry_price = starting_asset_price;
			// Sample the asset price at the expiry by sampling the process
			// along its random path, with small increments.
			for(unsigned int step = 0; step < n_steps; ++step)
				expiry_price *= std::exp((risk_free_interest_rate - 0.5 * std::pow(volatility, 2)) * time_step + volatility * rnd.Gauss(0, 1) * std::sqrt(time_step));
			call_profit += std::exp(-risk_free_interest_rate * expiry) * std::max(expiry_price - strike_price, 0.);
			put_profit  += std::exp(-risk_free_interest_rate * expiry) * std::max(strike_price - expiry_price, 0.);
		}
		call_profit /= block_size;
		call_estimates_sum += call_profit;
		call_squared_estimates_sum += std::pow(call_profit, 2);
		avg_call_option_price[i] = call_estimates_sum / (i + 1);
		put_profit /= block_size;
		put_estimates_sum += put_profit;
		put_squared_estimates_sum += std::pow(put_profit, 2);
		avg_put_option_price[i] = put_estimates_sum / (i + 1);
		if(i > 0)
		{
			std_call_option_price[i] = std::sqrt((call_squared_estimates_sum / (i + 1) - std::pow(avg_call_option_price[i], 2)) / i);
			std_put_option_price[i]  = std::sqrt((put_squared_estimates_sum / (i + 1) - std::pow(avg_put_option_price[i], 2)) / i);
		}
		else
		{
			std_call_option_price[i] = 0;
			std_put_option_price[i]  = 0;
		}
	}

	call_output_file.open("call-progressive-price.dat"),
	put_output_file.open("put-progressive-price.dat");
	for(unsigned int i = 0; i < n_blocks; ++i)
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
