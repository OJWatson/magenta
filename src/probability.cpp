#include "stdafx.h"
#include <random>
#include "probability.h"
#include <algorithm>
#include <cassert> // for error checking

using namespace std;

std::random_device rd;
//std::mt19937 generator (rd());
default_random_engine generator (rd());
std::uniform_real_distribution<double> uniform_0_1(0.0, 1.0);

//------------------------------------------------
// sample from uniform(a,b) distribution
double runif1(double a, double b) {
	std::uniform_real_distribution<double> uniform_a_b(a, b);
	return(uniform_a_b(generator));
}

//------------------------------------------------
// sample from uniform integer(a,b) distribution
int runiform_int_1(int a, int b) {
	std::uniform_int_distribution<int> unifrom_int_dist(a, b);
	return(unifrom_int_dist(generator));
}

//------------------------------------------------
// sample from given probability vector (that sums to pSum)
int sample1(std::vector<double> &p, double pSum) {
	double rand = pSum*uniform_0_1(generator);
	double z = 0;
	for (int i = 0; i<int(p.size()); i++) {
		z += p[i];
		if (rand<z)
			return i;
	}
	return(0);
}

//------------------------------------------------
// sample from gamma(alpha,beta) distribution
double rgamma1(double shape, double rate) {
	std::gamma_distribution<double> gamma_dist(shape, 1.0 / rate);
	double x = gamma_dist(generator);

	// check for zero or infinite values (corrects for bug in some compilers)
	while (x == 0 || (1.0 / x) == 0)
		x = gamma_dist(generator);

	return(x);
}

//------------------------------------------------
// draw from beta(alpha,beta) distribution
double rbeta1(double alpha, double beta) {
	double X1 = rgamma1(alpha, 1.0);
	double X2 = rgamma1(beta, 1.0);
	return(X1 / (X1 + X2));
}

//------------------------------------------------
// draw from univariate normal distribution
double rnorm1(double mean, double sd) {
	std::normal_distribution<double> normal_dist(mean, sd);
	return(normal_dist(generator));
}

//------------------------------------------------
// draw from lognormal distribution
double rlognorm1(double mean, double sd) {
	std::lognormal_distribution<double> lognormal_dist(mean, sd);
	return(lognormal_dist(generator));
}

//------------------------------------------------
// draw from poisson distribution
int rpoisson1(double mean) {
	std::poisson_distribution<int> poisson_dist(mean);
	return(poisson_dist(generator));
}

//------------------------------------------------
// draw from binomial distribution
int rbinomial1(int trials, double p) {
	std::binomial_distribution<int> binomial_dist(trials, p);
	return(binomial_dist(generator));
}

//------------------------------------------------
// draw from exponenetial distribution
double rexpdouble1(double lambda) {
	std::exponential_distribution<double> exponential_dist(lambda);
	return(exponential_dist(generator));
}

//------------------------------------------------
// draw from exponenetial distribution
int rexpint1(double lambda) {
	std::exponential_distribution<double> exponential_dist(lambda);
	return(static_cast<int>(exponential_dist(generator)));
}

//------------------------------------------------
// sample without replacement from vector
std::vector<int> sample_without_replacement(std::vector<int> &v, int c) {

	// catch incorrect vector vs sample size
	assert(static_cast<unsigned int>(c) < v.size() && "Sample size bigger than vector");

	// Create return vector
	std::vector<int> samples{ c };

	// Based on Knuth sample algorithm
	int n = c;
	int N = v.size();

	int t = 0; // total input records dealt with
	int m = 0; // number of items selected so far
	double u;

	while (m < n)
	{
		u = uniform_0_1(generator); // call a uniform(0,1) random number generator

		if ((N - t)*u >= n - m)
		{
			t++;
		}
		else
		{
			samples[m] = t;
			t++; m++;
		}
	}
	return(samples);
}

