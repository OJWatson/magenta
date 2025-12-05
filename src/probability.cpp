
#include "probability.h"

// --------------------------------------------------------------------------------------------------------------------------------
// Predefined distibutions with fixed means and one trial only
// --------------------------------------------------------------------------------------------------------------------------------

//------------------------------------------------
// draw from binomial distribution for cD
bool rbernoulli_cD() {
    return(R::unif_rand()<0.0676909);
} 

//------------------------------------------------
// draw from binomial distribution for cT
bool rbernoulli_cT() {
    return(R::unif_rand()<0.02179647);
}

//------------------------------------------------
// draw from binomial distribution for cU
bool rbernoulli_cU() {
    return(R::unif_rand()<0.006203);
}

// --------------------------------------------------------------------------------------------------------------------------------

//------------------------------------------------
// sample from uniform(0,1) distribution
double runif0_1() {
  return(R::runif(0, 1));
}

//------------------------------------------------
// sample from uniform(a,b) distribution
double runif1(double a, double b) {
    return(R::runif(a, b));
}

//------------------------------------------------
// sample from uniform integer(a,b) distribution
int runiform_int_1(int a, int b) {
    return(static_cast<int>(R::runif(a,b+1)));
}


//------------------------------------------------
// sample from given probability vector (that sums to pSum)
int sample1(const std::vector<double> &p, double pSum) {
    double rand = pSum*R::unif_rand();
    double z = 0;
    for (int i = 0; i<int(p.size()); i++) {
        z += p[i];
        if (rand < z)
            return i;
    }
    return(0);
}

//------------------------------------------------
// sample from given probability vector (that sums to pSum)
int sample1_ints(const std::vector<int> &p, int pSum) {
    double rand = R::unif_rand()*pSum;
    double z = 0;
    for (int i = 0; i<int(p.size()); i++) {
        z += p[i];
        if (rand < z)
            return i;
    }
    return(0);
}



//------------------------------------------------
// sample given sorted random vector, cumsum vector, and n to sammple
void samplerandoms(std::vector<double> &r, std::vector<double> &p, int n, std::vector<int> &bite_storage_queue) {
  
  int n_count = 0;
  int iteration = 0;
  while (n_count < n) 
  {
    if (r[n_count] <= p[iteration]) {
      bite_storage_queue.emplace_back(iteration);
      n_count++;
      
      while (r[n_count]<p[iteration]) {
        bite_storage_queue.emplace_back(iteration);
        n_count++;
      }
    }
    iteration++;
  }
  
}

//------------------------------------------------
// sample from gamma(alpha,beta) distribution
double rgamma1(double shape, double rate) {

    double x = R::rgamma(shape, 1.0 / rate);
    // check for zero or infinite values (corrects for bug in some compilers)
    while (x == 0 || (1.0 / x) == 0)
    {
        x = R::rgamma(shape, 1.0 / rate);
    }
    
    return(x);
}

//------------------------------------------------
// sample from nbinom(size,prob) distribution
double rnbinom_mu(double size, double mu)
{
  return (mu == 0) ? 0 : rpoisson1(rgamma1(size, mu / size));
}

//------------------------------------------------
// draw from beta(alpha,beta) distribution
double rbeta1(double alpha, double beta) {
    return(R::rbeta(alpha, beta));
}

//------------------------------------------------
// draw from univariate normal distribution
double rnorm1(double mean, double sd) {
    return(R::rnorm(mean,sd));
}

//------------------------------------------------
// draw from lognormal distribution
double rlognorm1(double mean, double sd) {
    return(R::rlnorm(mean,sd));
}

//------------------------------------------------
// draw from poisson distribution
int rpoisson1(double mean) {
    return(R::rpois(mean));
}

//------------------------------------------------
// draw from binomial distribution
int rbinomial1(int trials, double p) {
    return(R::rbinom(trials,p));
}

// Simple implementation of the rmultinom.c that emplaces the draws to the vector provided
// Optimized: replaced individual emplace_back loops with batch insert operations
void rmultinomN(int n, std::vector<double> &prob, double p_tot, int K, std::vector<int> &output)
{
  int k;
  int draw;
  
  for (k = 0; k < K - 1; k++) { /* (p_tot, n) are for "remaining binomial" */
    
    draw = rbinomial1(n, (prob[k] / p_tot));
    n -= draw;
    
    // Batch insert instead of individual emplace_back calls
    if (draw > 0) {
      output.insert(output.end(), draw, k);
    }
    
    if (n <= 0) /* we have all*/ return;
    p_tot -= prob[k]; /* i.e. = sum(prob[(k+1):K]) */
  }
  
  // Batch insert remaining draws
  if (n > 0) {
    output.insert(output.end(), n, K - 1);
  }
  return;
}

//------------------------------------------------
// draw from smarter binomial distribution
int rbinomialsmarter1(unsigned int trials, double p) {
  
  double q = pow(1 - p, trials);
  if (rbernoulli1(q)) {
    return(0);
  }
  double r = 1 - q; // remaining mass
  double tmp1 = p / (1 - p);
  for (unsigned int i = 0; i<trials; i++) {
    q *= (trials - i) / double(i + 1) * tmp1;
    if (rbernoulli1(q / r)) {
      return(i + 1);
    }
    r -= q;
  }
  return(-1);
  
}

//------------------------------------------------
// draw from bernoulli distribution
bool rbernoulli1(double p) {
    return(R::unif_rand()<p);
}

//------------------------------------------------
// draw from exponenetial distribution
double rexpdouble1(double lambda) {
    return(R::rexp(lambda));
}

//------------------------------------------------
// draw from exponenetial distribution
int rexpint1(double lambda) {
    return(static_cast<int>(R::rexp(lambda)));
}

//------------------------------------------------
// sample without replacement from vector
std::vector<int> sample_without_replacement(std::vector<int> &v, int c) {
    
    if(static_cast<unsigned int>(c) > v.size()) {
      return(v);
    } else {
    
    // Create return vector
    std::vector<int> samples( c );
    
    // Based on Knuth sample algorithm
    int n = c;
    int N = v.size();
    
    int t = 0; // total input records dealt with
    int m = 0; // number of items selected so far
    double u;
    
    while (m < n)
    {
        u = R::unif_rand(); // call a uniform(0,1) random number generator
        
        if ((N - t)*u >= n - m)
        {
            t++;
        }
        else
        {
            samples[m] = v[t];
            t++; m++;
        }
    }
    return(samples);
    }
}

//------------------------------------------------
// sample with replacement from vector
std::vector<int> sample_with_replacement(std::vector<int> &v, int c) {
    
    // Create return vector
    std::vector<int> samples( c );
    int N = v.size();
    
    for(int i=0 ; i<c ; i++)
    { 
        samples[i] = v[runiform_int_1(0, N)]; // call a uniform(0,1) random integer
    }
    
    return(samples);
}

//------------------------------------------------
// sample with replacement from max integer
std::vector<int> sample_with_replacement_from_max(int max, int c) {
    
    // Create return vector
    std::vector<int> samples(c);
    
    for (int i=0; i<c; i++)
    {
        samples[i] = runiform_int_1(0, max); // call a uniform(0,N) random intenger
    }
    
    return(samples);
}

//------------------------------------------------
// shuffle integer vector
inline int randWrapper(const int n) { return std::floor(R::unif_rand()*n); }

void shuffle_integer_vector(std::vector<int> &vec) {
  
  int vec_size = vec.size();
  for (int i = vec_size-1; i > 0; i--) {
    std::swap(vec[i], vec[std::floor(R::unif_rand()*i)]);
  }

}


//------------------------------------------------
// hill function
double hill_function(double x, double n, double kA) {
  return((pow(x, n)) / (pow(x, n) + pow(kA, n)));
}


//------------------------------------------------
// is integer zero
bool is_int_zero(int i) {
  return i == 0;
}


//------------------------------------------------
// Is double zero function
bool is_double_zero(double i) {
  return i == 0.0;
};
