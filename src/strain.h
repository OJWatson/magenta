//
//  magenta
//  strain.h
//
//  Created: OJ on 13/01/2017
//	Most recent edits: OJ on 10/03/2017
//
//  Distributed under the MIT software licence - see Notes.c file for details
//
//  Class describing a parasite strain.
//
// ---------------------------------------------------------------------------


#ifndef STRAIN_H
#define STRAIN_H
#include <vector>
#include <queue>
#include <bitset>
#include <random>
#include "probability.h"
#include "parameters.h"
#include "util.h"
#include <boost/utility.hpp>

class Strain {
  
public:
  enum InfectionStatus
  {
    SUSCEPTIBLE,	// 0
    DISEASED,		// 1
    ASYMPTOMATIC,	// 2
    SUBPATENT,		// 3
    TREATED,		// 4
    PROPHYLAXIS,	// 5
    NUMBER_OF_STATES = 6
  };

  // Infection transition options. 
  const static std::vector<InfectionStatus> m_transition_vector;
  
  // temporary barcodes for all purposes
  static boost::dynamic_bitset<> temp_barcode;
  static std::vector<boost::dynamic_bitset<> > temp_barcode_pair;
  static boost::dynamic_bitset<> temp_identity_barcode;
  static boost::dynamic_bitset<> temp_crossovers;
  
private:
  
  boost::dynamic_bitset<> m_barcode;									// barcode sequence
  InfectionStatus m_strain_infection_status;					// infection status associated with a strain
  int m_day_of_strain_infection_status_change;				// day that strain would move infection status
  int m_day_of_acquisition;                           // day that the strain was acquired by an individual/mosquito

public:
  
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // CONSTRUCTORS & DESTRUCTOR
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  // Default class constructor which will inialise a random strain
  Strain();
  
  // Known barcode constructor
  Strain(boost::dynamic_bitset<> barcode);
  
  // Known barcode and infection status constructor
  Strain(boost::dynamic_bitset<> barcode, const Strain::InfectionStatus &infectionStatus);
  
  // Known barcode, infection status and day of infection status change outcome constructor
  Strain(boost::dynamic_bitset<> barcode, const Strain::InfectionStatus &infectionStatus, const int &dayOfInfectionStatusChange);
  
  // Known barcode, infection status and day of infection status change outcome constructor and day of acquisition
  Strain(boost::dynamic_bitset<> barcode, const Strain::InfectionStatus &infectionStatus, 
         const int &dayOfInfectionStatusChange, const int & day_of_acquisition);
  
  // Destructor
  // ~Strain();
  
  
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // GETTERS
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  // Get barcode
  boost::dynamic_bitset<> get_m_barcode() const { return(m_barcode); }
  
  // Get strain infection status
  InfectionStatus get_m_strain_infection_status() { return(m_strain_infection_status); }
  
  // Get day of strain infection status
  int get_m_day_of_strain_infection_status_change() { return(m_day_of_strain_infection_status_change); }
  
  // Get day of strain acquisition
  int get_m_day_of_strain_acquisition() { return(m_day_of_acquisition); }
  
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // MEMBER FUNCTIONS
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  // Work out the strain's relative onward contribution given the fitness costs in the parameters object
  double relative_contribution(const Parameters &parameters);
  
  // Work out if the strain caused late parasitological failure
  double late_paristological_failure_prob(const Parameters &parameters, int drug_choice);
  
  // Work out if the strain is resistant at all for the drug given
  bool resistant_boolean(const Parameters &parameters, int drug_choice);
  
  // Work out if the strain is vector adapted
  bool vector_adapted_boolean(const Parameters &parameters);

  // Work out if the strain is resitsant at any resistant loci
  bool resistant_at_any_loci_boolean(const Parameters &parameters);
  
  // Get barcode position
  bool barcode_position(unsigned int position);
  
  // to_ulong but for positions in a bitset
  unsigned long bitset_positions_to_ulong(const std::vector<unsigned int> &bit_positions);
  
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // SETTERS
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  // Set barcode pointer
  void set_m_barcode(boost::dynamic_bitset<> x) { m_barcode = x; }
  
  // Set strain infection status
  void set_m_strain_infection_status(InfectionStatus x) { m_strain_infection_status = x; }
  
  // Set strain infection status day change
  void set_m_day_of_strain_infection_status_change(int x) { m_day_of_strain_infection_status_change = x; }
  
  // Set strain acquisition day
  void set_m_day_of_strain_acquisition(int x) { m_day_of_acquisition = x; }
  
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // CLASS STATIC FUNCTIONS
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  // Generate next barcode. For non IBD this is just a random. IBD its the next identity
  static boost::dynamic_bitset<> generate_next_barcode();
  
  // Generate a random barcode
  static boost::dynamic_bitset<> generate_random_barcode();
  
  // Generate a random ordinary barcode, i.e. 1 bit per loci 
  static boost::dynamic_bitset<> generate_random_ordinary_barcode();
  
  // Generate a random ordinary barcode, i.e. 1 bit per loci 
  static boost::dynamic_bitset<> generate_random_ibd_barcode();
  
  // Generate a random identity barcode, i.e. where the barcode represents num_loci * ibd_length
  static boost::dynamic_bitset<> generate_random_identity_barcode();
  
  // Generate the next identity barcode, i.e. if it was 0101 it will now be 1111, if num_loci = 2
  static boost::dynamic_bitset<> generate_next_ibd_barcode();
  
  // Generate a random recombinant barcode given two barcodes
  static std::vector<boost::dynamic_bitset<> > generate_recombinant_barcode(boost::dynamic_bitset<> x, boost::dynamic_bitset<> y);
  
  // Generate a random barcode given probability of each SNP, i.e. PLAF
  static boost::dynamic_bitset<> generate_random_barcode_given_SNP_frequencies(std::vector<double> &x);
  
  // Generate a dependent barcode given probability of each SNP, i.e. PLAF
  static boost::dynamic_bitset<> generate_dependent_barcode_given_SNP_frequencies(std::vector<double> &x);
  
  // Create a bitset by replicating each bit n times
  static boost::dynamic_bitset<> replicate_by_bit(boost::dynamic_bitset<> x, unsigned int n);
  
  // Turns our ibd barcode into a vector of the ints making it up
  static std::vector<unsigned int> ibd_barcode_to_integer_vector(boost::dynamic_bitset<> x);
  
  // distances between one bitset and a vector range 
  static unsigned int distance_of_bitset_a_and_x(boost::dynamic_bitset<> a, 
                                                  std::vector<boost::dynamic_bitset<> >::const_iterator start, 
                                                  std::vector<boost::dynamic_bitset<> >::const_iterator end);
  
  // distances between one bitset and a vector range of bitsets
  static unsigned int ibd_distance_of_bitset_a_and_x(boost::dynamic_bitset<> a, 
                                                      std::vector<boost::dynamic_bitset<> >::const_iterator start, 
                                                      std::vector<boost::dynamic_bitset<> >::const_iterator end);
  
  // distances between one bitset and a vector of vectors of bitsets
  static double distance_of_bitset_a_and_vec_x(boost::dynamic_bitset<> a, 
                                                std::vector<std::vector< boost::dynamic_bitset<> > >::const_iterator start, 
                                                std::vector<std::vector< boost::dynamic_bitset<> > >::const_iterator end);
  
  // mean distance between all bitsets within a vector of bitsets
  static double distance_mean_within_bitsets(std::vector<boost::dynamic_bitset<> > x, unsigned int bl, bool unique_only = false);

  // mean ibd distance between all bitsets within a vector of bitsets
  static double ibd_distance_mean_within_bitsets(std::vector<boost::dynamic_bitset<> > x, bool unique_only = false);
  
  // OTHER STATICS
  
  // barcode true at all locations specified by vector<int>
  static bool all_at_positions(boost::dynamic_bitset<> x, std::vector<unsigned int> pos);
  
  // barcode true at any locations specified by vector<int>
  static bool any_at_positions(boost::dynamic_bitset<> x, std::vector<unsigned int> pos);
  
};

#endif
