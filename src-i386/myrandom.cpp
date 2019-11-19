// simple wrappers around pseudo random number generator in C++11
// and normal distribution

#include <random>
#include "myrandom.h"

# define M_PI           3.14159265358979323846  /* pi */

// set up random number generator
std::random_device rd;  //Will be used to obtain a seed for the random number engine
std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()


double runif()
{
// Uniform random number generator, see https://en.cppreference.com/w/cpp/numeric/random/uniform_real_distribution
	std::uniform_real_distribution<double> dis(0, 1);
	return dis(gen);
}

double rnorm()
{
	// N(0,1) random number generator
	std::normal_distribution<> d{0,1};
	return d(gen);
}

double dnorm(double x, double mean, double sigma)
{
	return 1/(sigma*sqrt(2*M_PI)) * exp( -pow((x-mean)/sigma, 2) / 2 );
}

double dnorm_std(double x)
{
	return 1/sqrt(2*M_PI) * exp( -x*x/2 );
}

double rexp(double rate)
{
	std::exponential_distribution<> d(rate);
	return d(gen);
}