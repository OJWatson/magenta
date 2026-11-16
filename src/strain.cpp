#include "strain.h"

// ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// CONSTRUCTORS
// ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

// Default class constructor which will inialise a random strain
Strain::Strain() :
  m_barcode(Strain::generate_next_barcode())
{
}

// Known barcode constructor
Strain::Strain(boost::dynamic_bitset<> barcode) :
  m_barcode(barcode)
{
}

// Known barcode and infection status constructor
Strain::Strain(boost::dynamic_bitset<> barcode, const Strain::InfectionStatus &infectionStatus) :
  m_barcode(barcode),
  m_strain_infection_status(infectionStatus)
{
}

// Known barcode, infection status and day of infection status change outcome constructor
Strain::Strain(boost::dynamic_bitset<> barcode, const Strain::InfectionStatus &infectionStatus, const int &dayOfInfectionStatusChange) :
  m_barcode(barcode),
  m_strain_infection_status(infectionStatus),
  m_day_of_strain_infection_status_change(dayOfInfectionStatusChange)
{
}

// Known barcode, infection status and day of infection status change outcome constructor and acquisition dat
Strain::Strain(boost::dynamic_bitset<> barcode, const Strain::InfectionStatus &infectionStatus, 
               const int &dayOfInfectionStatusChange, const int & day_of_acquisition) :
  m_barcode(barcode),
  m_strain_infection_status(infectionStatus),
  m_day_of_strain_infection_status_change(dayOfInfectionStatusChange),
  m_day_of_acquisition(day_of_acquisition)
{
}

// ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// DESTRUCTOR
// ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------



// ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// STATIC INITIALIZERS
// ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

// Initialise static class const vectors
const std::vector<Strain::InfectionStatus> Strain::m_transition_vector{ SUSCEPTIBLE, DISEASED, ASYMPTOMATIC, SUBPATENT, TREATED, PROPHYLAXIS };

// temp barcodes for speed
boost::dynamic_bitset<> Strain::temp_barcode(Parameters::g_barcode_length);
std::vector<boost::dynamic_bitset<> > Strain::temp_barcode_pair(2);
boost::dynamic_bitset<> Strain::temp_identity_barcode(Parameters::g_ibd_length);
boost::dynamic_bitset<> Strain::temp_crossovers(Parameters::g_num_loci);

// ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// MEMBER FUNCTIONS
// ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

// Work out the strain's relative onward contribution given the fitness costs in the parameters object
double Strain::relative_contribution(const Parameters &parameters){
  // the cost table is indexed by the bits at the resistance loci, bit j being
  // the j-th locus, which is how drug_list_create() builds it
  return(parameters.g_cost_of_resistance[bitset_positions_to_ulong(parameters.g_resistance_loci)]);
}

// Work out if the strain caused late parasitological failure
double Strain::late_paristological_failure_prob(const Parameters &parameters, int drug_choice){
  
  return(1 - parameters.g_drugs[drug_choice].get_prob_of_lpf_barcode(m_barcode));
  
}

// Work out if the strain is resistant at all for the drug given
bool Strain::resistant_boolean(const Parameters &parameters, int drug_choice){
  
  if(parameters.g_drugs[drug_choice].get_prob_of_lpf_barcode(m_barcode) == 1.0) { 
    return(false);
  } else {
    return(true);
  }
}

// Get barcode position
bool Strain::barcode_position(unsigned int position){
  return(m_barcode[position]);
}

// to_ulong but for positions in a bitset
unsigned long Strain::bitset_positions_to_ulong(const std::vector<unsigned int> &bit_positions){
  
  unsigned long mask = 1;
  unsigned long result = 0;
  for (unsigned int start_bit = 0 ; start_bit < bit_positions.size(); start_bit++) {
    if (m_barcode[bit_positions[start_bit]]) {
      result |= mask;
    }
    mask <<= 1;
  }
  return result;
  
}

// Work out if the strain is vector adapted
bool Strain::vector_adapted_boolean(const Parameters &parameters){
  
  bool adapted = true;
  for(auto v : parameters.g_vector_adaptation_loci) {
    if(!m_barcode[v]){
      adapted = false;
    }
  }
  return(adapted);
  
}

// Work out if the strain is resitsance at any resistant loci
bool Strain::resistant_at_any_loci_boolean(const Parameters &parameters){
  
  for(auto v : parameters.g_resistance_loci) {
    if(m_barcode[v]){
      return(true);
    }
  }
  return(false);
  
}




// ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// STRAIN - UTILS
// ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

// Generate next barcode. For non IBD this is just a random. IBD its the next identity
boost::dynamic_bitset<> Strain::generate_next_barcode()
{
  // Match barcode type
  switch (Parameters::g_barcode_type)
  {
  case Parameters::ORDINARY:
    // # is it drawn non independently from the island
    if(Parameters::g_island_imports_plaf_linked_flag) {
      return(generate_dependent_barcode_given_SNP_frequencies(Parameters::g_plaf));  
    } else {
      return(generate_random_barcode_given_SNP_frequencies(Parameters::g_plaf));
    }
  case Parameters::IBD:
    return(generate_next_ibd_barcode());
  default:
    Rcpp::stop("Unrecognised barcode_type");
  break;
  }
  
}

// Generate a random barcode fiven plaf
boost::dynamic_bitset<> Strain::generate_random_barcode()
{
  // Match barcode type
  switch (Parameters::g_barcode_type)
  {
  case Parameters::ORDINARY:
    return(generate_random_barcode_given_SNP_frequencies(Parameters::g_plaf));
  case Parameters::IBD:
    return(generate_random_ibd_barcode());
  default:
    Rcpp::stop("Unrecognised barcode_type");
  break;
  }
}

// Generate a random ordinary barcode, i.e. 1 bit per loci 
boost::dynamic_bitset<> Strain::generate_random_ordinary_barcode()
{
  for(unsigned int i = 0; i < Parameters::g_barcode_length ; i++){
    Strain::temp_barcode[i] = rbernoulli1(0.5);
  }
  return(Strain::temp_barcode);
}

// Generate a random identity barcode, i.e. bitset of length equal to ibd_length
boost::dynamic_bitset<> Strain::generate_random_identity_barcode()
{
  for(unsigned int i = 0; i < Parameters::g_ibd_length ; i++){
    Strain::temp_identity_barcode[i] = rbernoulli1(0.5);
  }
  
  return(Strain::temp_identity_barcode);
}

// Generate a random identity barcode, i.e. where the barcode represents num_loci * ibd_length
boost::dynamic_bitset<> Strain::generate_random_ibd_barcode()
{
  unsigned int count = 0;
 
  // create and fill the temp barcode
  for(unsigned int i = 0; i < Parameters::g_num_loci ; i++) {
    
    Strain::temp_identity_barcode = generate_random_identity_barcode();
    
    for(unsigned int j = 0; j < Parameters::g_ibd_length ; j++) {
      
      Strain::temp_barcode[count++] =  Strain::temp_identity_barcode[j];
    }
  }
  
  return(Strain::temp_barcode);
}

// Generate a random identity barcode, i.e. where the barcode represents num_loci * ibd_length
boost::dynamic_bitset<> Strain::generate_next_ibd_barcode()
{
  // create the next identity bitset
  Strain::temp_identity_barcode = boost::dynamic_bitset<>(Parameters::g_ibd_length, Parameters::g_identity_id++);
  
  unsigned int count = 0;
  
  // then fill in th temp by repeating this
  for(unsigned int i = 0; i < Parameters::g_num_loci ; i++) {
    
    for(unsigned int j = 0; j < Parameters::g_ibd_length ; j++) {
      Strain::temp_barcode[count++] =  Strain::temp_identity_barcode[j];
    }
  }
  
  return(Strain::temp_barcode);
}

// Generate a crossovers bitset using the prob_crossovers, and is thus length of num_loci
boost::dynamic_bitset<> generate_crossover()
{
  for(unsigned int i = 0; i < Parameters::g_num_loci ; i++){
    Strain::temp_crossovers[i] = rbernoulli1(Parameters::g_prob_crossover[i]);
  }
  return(Strain::temp_crossovers);
}

// Generate a random recombinant barcode given two barcodes
std::vector<boost::dynamic_bitset<> > Strain::generate_recombinant_barcode(boost::dynamic_bitset<> x, boost::dynamic_bitset<> y)
{
  // find the different positions, then where these cross with another 
  // bitset created using the prob_crossover which represents chance of segragation. 
  // i.e. 0.5 equals indepedent segregation, 0 equals loci next to each other and
  // 0% recombination. Then where these differ with x. 
  // For the IBD style, our random barcode
  
  // When you have done the switch then set the first bit of ibd_block length equal to the parent's
  // first bit
  
  // Match barcode type
  switch (Parameters::g_barcode_type)
  {
  case Parameters::ORDINARY:
    
    Strain::temp_barcode_pair[0] = (((x ^ y) & generate_random_barcode_given_SNP_frequencies(Parameters::g_prob_crossover) ) ^ x);
    Strain::temp_barcode_pair[0][0] = x[0];
    Strain::temp_barcode_pair[1] = (((x ^ y) & generate_random_barcode_given_SNP_frequencies(Parameters::g_prob_crossover) ) ^ y);
    Strain::temp_barcode_pair[1][0] = y[0];
    
    return(Strain::temp_barcode_pair);
    
  case Parameters::IBD:
    
    Strain::temp_barcode_pair[0] = (((x ^ y) & replicate_by_bit(generate_crossover(),Parameters::g_ibd_length) ) ^ x);
    for(unsigned int j = 0; j < Parameters::g_ibd_length ; j++) {
      Strain::temp_barcode_pair[0][j] = x[j];
    }
    
    Strain::temp_barcode_pair[1] = (((x ^ y) & replicate_by_bit(generate_crossover(),Parameters::g_ibd_length) ) ^ y);
    for(unsigned int j = 0; j < Parameters::g_ibd_length ; j++) {
      Strain::temp_barcode_pair[1][j] = y[j];
    }
    
    return(Strain::temp_barcode_pair);
    
  default:
    Rcpp::stop("Unrecognised barcode_type");
  break;
  }

}

// Generate a random barcode given probability of each SNP, i.e. PLAF
boost::dynamic_bitset<> Strain::generate_random_barcode_given_SNP_frequencies(std::vector<double> &x)
{
  
  // create and fill the temp barcode
  for(unsigned int i = 0; i < Parameters::g_barcode_length ; i++){
    Strain::temp_barcode[i] = rbernoulli1(x[i]);
  }

  // return barcode
  return(Strain::temp_barcode);
}

// Generate a random barcode given probability of each SNP, i.e. PLAF
boost::dynamic_bitset<> Strain::generate_dependent_barcode_given_SNP_frequencies(std::vector<double> &x)
{
  
    auto it = find_if_not(x.begin(), x.end(), is_double_zero);
    int first_nonzero_idx = distance(x.begin(), it);

    if(rbernoulli1(x[first_nonzero_idx])) {
      
      // create and fill the temp barcode with plaf positions not zero as 1
      for(unsigned int i = 0; i < Parameters::g_barcode_length ; i++){
        Strain::temp_barcode[i] = !is_double_zero(x[i]);
      }
      
    } else {
      
      // create and fill the temp barcode with 0
      for(unsigned int i = 0; i < Parameters::g_barcode_length ; i++){
        Strain::temp_barcode[i] = false;
      }
      
    }
    
  // return barcode
  return(Strain::temp_barcode);
  
}


// Stretches a bitset by replicating each bit n times
boost::dynamic_bitset<> Strain::replicate_by_bit(boost::dynamic_bitset<> x, unsigned int n)
{
  
  unsigned int count = 0;
  
  // create and fill the temp barcode
  for(unsigned int i = 0; i < Parameters::g_num_loci ; i++) {
    
    for(unsigned int j = 0; j < n ; j++) {
      
      Strain::temp_barcode[count++] = x[i];
      
    }
  }
  
  // return barcode
  return(Strain::temp_barcode);
}

// Turns our ibd barcode into a vector of the ints making it up
std::vector<unsigned int> Strain::ibd_barcode_to_integer_vector(boost::dynamic_bitset<> x)
{
  
  unsigned int mask = 1;
  unsigned int result = 0;
  std::vector<unsigned int> vec(Parameters::g_num_loci);
  
  for(unsigned int j = 0; j < Parameters::g_num_loci; j++)
  {
    mask = 1;
    result = 0;
  for (size_t i = (0+((j)*Parameters::g_ibd_length)); i < ((j+1)*Parameters::g_ibd_length); ++ i) {
    if (x.test(i))
      result |= mask;
    mask <<= 1;
  }
  vec[j] = result;
  }
  // return barcode
  return(vec);
}

// --------------------------------------------------------------------------------------------------
// DISTANCES
// --------------------------------------------------------------------------------------------------

// distances between one bitset and a vector range of bitsets
unsigned int Strain::distance_of_bitset_a_and_x(boost::dynamic_bitset<> a, 
                                 std::vector<boost::dynamic_bitset<> >::const_iterator start, 
                                 std::vector<boost::dynamic_bitset<> >::const_iterator end)
{
  
  unsigned int distance = 0;
  while (start != end) // while it hasn't reach the end
  {
    distance = distance + ((~(a ^ *start)).count());
    ++start; // and iterate to the next element
  }
  // return barcode
  return(distance);
}

// distances between one bitset and a vector range of bitsets
unsigned int Strain::ibd_distance_of_bitset_a_and_x(boost::dynamic_bitset<> a, 
                                                std::vector<boost::dynamic_bitset<> >::const_iterator start, 
                                                std::vector<boost::dynamic_bitset<> >::const_iterator end)
{
  
  unsigned int distance = 0;
  std::vector<unsigned int> vec(Parameters::g_num_loci);
  
  while (start != end) // while it hasn't reach the end
  {
    vec = Strain::ibd_barcode_to_integer_vector(a ^ *start);
    for(auto c : vec){
      if(!c) distance++;
    }
    ++start; // and iterate to the next element
  }
  // return barcode
  return(distance);
}

// distances between one bitset and a vector of vectors of bitsets
double Strain::distance_of_bitset_a_and_vec_x(boost::dynamic_bitset<> a, 
                                                std::vector<std::vector< boost::dynamic_bitset<> > >::const_iterator start, 
                                                std::vector<std::vector< boost::dynamic_bitset<> > >::const_iterator end)
{
  
  double distance = 0.0;
  while (start != end) // while it hasn't reach the end
  {
    distance = distance + (
      (Strain::distance_of_bitset_a_and_x(a, (*start).begin(), (*start).end())) / 
      ( (*start).end() - (*start).begin() )
    );
    ++start; // and iterate to the next element
  }
  // return barcode
  return(distance);
}

// mean distance between all bitsets within a vector of bitsets
double Strain::distance_mean_within_bitsets(std::vector<boost::dynamic_bitset<> > x, unsigned int bl, bool unique_only)
{
  
  unsigned int distance = 0;
  
  // should we unique the bitsets first
  if(unique_only) {
    std::sort( x.begin(), x.end() );
    x.erase( unique( x.begin(), x.end() ), x.end() );
  }
  
  for(unsigned int i = 1; i < (x.size()); i++)
  {
    distance = distance + Strain::distance_of_bitset_a_and_x(x[i-1], x.begin()+i, x.end());
  }
  
  // return distance mean, i.e. divided by the length of the barcode and the nC2 combinations that were calculated
  return(distance / (bl * (x.size() * (x.size() - 1) * 0.5) ) );
}

// mean distance between all bitsets within a vector of bitsets
double Strain::ibd_distance_mean_within_bitsets(std::vector<boost::dynamic_bitset<> > x, bool unique_only)
{
  
  unsigned int distance = 0;
  

  // should we unique the bitsets first
  if(unique_only) {
    std::sort( x.begin(), x.end() );
    x.erase( unique( x.begin(), x.end() ), x.end() );
  }
  
  for(unsigned int i = 1; i < (x.size()); i++)
  {
    distance = distance + Strain::ibd_distance_of_bitset_a_and_x(x[i-1], x.begin()+i, x.end());
  }
  
  // return distance mean, i.e. divided by the num of loci and the nC2 combinations that were calculated
  return(distance / (Parameters::g_num_loci * (x.size() * (x.size() - 1) * 0.5) ) );
}

// --------------------------------------------------------------------------------------------------
// MISC
// --------------------------------------------------------------------------------------------------


// barcode true at all locations specified by vector<int>
bool Strain::all_at_positions(boost::dynamic_bitset<> x, std::vector<unsigned int> pos)
{
  bool all = true;
  for(auto v : pos) {
    if(!x[v]){
      all = false;
    }
  }
  return(all);
  
}

// barcode true at any locations specified by vector<int>
bool Strain::any_at_positions(boost::dynamic_bitset<> x, std::vector<unsigned int> pos)
{
  for(auto v : pos) {
    if(x[v]){
      return(true);
    }
  }
  return(false);
  
}

// --------------------------------------------------------------------------------------------------
// TESTS
// --------------------------------------------------------------------------------------------------

// PLAF test
// [[Rcpp::export]]
SEXP test_barcode_from_PLAF(Rcpp::NumericVector plaf, unsigned int n)
{
  Parameters::g_barcode_length = n;
  Strain::temp_barcode = boost::dynamic_bitset<>(Parameters::g_barcode_length);
  std::vector<double> plaf_c = Rcpp::as<std::vector<double> >(plaf);
  boost::dynamic_bitset<> b = Strain::generate_random_barcode_given_SNP_frequencies(plaf_c);
  return(bitset_to_sexp(b, n));
}

// PLAF dependent test
// [[Rcpp::export]]
SEXP test_dependent_barcode_from_PLAF(Rcpp::NumericVector plaf, unsigned int n)
{
  Parameters::g_barcode_length = n;
  Strain::temp_barcode = boost::dynamic_bitset<>(Parameters::g_barcode_length);
  std::vector<double> plaf_c = Rcpp::as<std::vector<double> >(plaf);
  boost::dynamic_bitset<> b = Strain::generate_dependent_barcode_given_SNP_frequencies(plaf_c);
  return(bitset_to_sexp(b, n));
}

// IBD recombinant test
// [[Rcpp::export]]
SEXP test_recombinant_with_ibd(SEXP barcode_1,
                               SEXP barcode_2,
                               unsigned int bl, unsigned int nl,
                               unsigned int ib, Rcpp::NumericVector pc
                               )
{
  Parameters::g_barcode_length = bl;
  Parameters::g_num_loci = nl;
  Parameters::g_ibd_length = ib;
  Parameters::g_prob_crossover = Rcpp::as<std::vector<double> >(pc);
  Parameters::g_barcode_type = Parameters::IBD;
  
  Strain::temp_barcode = boost::dynamic_bitset<>(Parameters::g_barcode_length);
  Strain::temp_crossovers = boost::dynamic_bitset<>(Parameters::g_num_loci);
  
  boost::dynamic_bitset<> barcode_a = sexp_to_bitset(barcode_1, bl);
  boost::dynamic_bitset<> barcode_b = sexp_to_bitset(barcode_2, bl);
  std::vector<boost::dynamic_bitset<> > temp_barcodes = Strain::generate_recombinant_barcode(barcode_a, barcode_b);
  
  return(bitset_to_sexp(temp_barcodes[1], bl));
}

// generate next ibd barcode
// [[Rcpp::export]]
SEXP test_generate_next_ibd(unsigned int bl, unsigned int nl,
                               unsigned int ib, Rcpp::NumericVector pc,
                               unsigned int id
)
{
  Parameters::g_barcode_length = bl;
  Parameters::g_num_loci = nl;
  Parameters::g_ibd_length = ib;
  Parameters::g_prob_crossover = Rcpp::as<std::vector<double> >(pc);
  Parameters::g_barcode_type = Parameters::IBD;
  Parameters::g_identity_id = id;
  
  Strain::temp_barcode = boost::dynamic_bitset<>(Parameters::g_barcode_length);
  Strain::temp_crossovers = boost::dynamic_bitset<>(Parameters::g_num_loci);
  Strain::temp_identity_barcode = boost::dynamic_bitset<>(Parameters::g_ibd_length);
  boost::dynamic_bitset<> barcode_a = Strain::generate_next_ibd_barcode();
  return(bitset_to_sexp(barcode_a, bl));
}

// ibd to vector of ints
// [[Rcpp::export]]
Rcpp::List test_ibd_conversion(SEXP barcode, unsigned int bl, 
                        unsigned int nl,
                        unsigned int ib
)
{
  Parameters::g_barcode_length = bl;
  Parameters::g_num_loci = nl;
  Parameters::g_ibd_length = ib;
  Parameters::g_barcode_type = Parameters::IBD;
  
  Strain::temp_barcode = boost::dynamic_bitset<>(Parameters::g_barcode_length);
  
  boost::dynamic_bitset<> barcode_a = sexp_to_bitset(barcode, bl);
  std::vector<unsigned int> res(Strain::ibd_barcode_to_integer_vector(barcode_a));
  return(Rcpp::List::create(Rcpp::Named("vec") = res));
}


// fitness cost lookup test
// [[Rcpp::export]]
double test_relative_contribution(SEXP barcode, unsigned int bl,
                                  Rcpp::IntegerVector resistance_loci,
                                  Rcpp::NumericVector cost_of_resistance)
{
  Parameters::g_barcode_length = bl;
  Parameters parameters;
  parameters.g_resistance_loci = Rcpp::as<std::vector<unsigned int> >(resistance_loci);
  parameters.g_cost_of_resistance = Rcpp::as<std::vector<double> >(cost_of_resistance);
  Strain strain(sexp_to_bitset(barcode, bl));
  return(strain.relative_contribution(parameters));
}
