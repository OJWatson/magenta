//
//  magenta
//  parameters.h
//
//  Created: OJ on 17/01/2017
//
//  Distributed under the MIT software licence - see Notes.c file for details
//
//  Class describing all global parameters.
//
// ---------------------------------------------------------------------------

#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <cmath>
#include <vector>
#include <string>
#include <queue>
#include <bitset>
#include <boost/dynamic_bitset.hpp>
#include "drug.h"

class Parameters {
  
public:
  
  int g_current_time;
  int g_calendar_day;
  std::vector<double> g_theta;
  double g_years;
  
  // spatial
  enum g_spatial_type_enum{
    NON, //0
    ISLAND, //1
    METAPOPULATION, //2
    NUMBER_OF_SPATIAL_TYPE_OPTIONS = 3
  };
  
  // spatially relevant variables
  // metapopulation related
  Parameters::g_spatial_type_enum g_spatial_type;
  std::vector<int> g_spatial_imported_cotransmission_frequencies;
  std::vector<int> g_spatial_imported_oocyst_frequencies;
  std::vector<int> g_spatial_exported_cotransmission_frequencies;
  std::vector<int> g_spatial_exported_oocyst_frequencies;
  unsigned int g_spatial_total_exported_barcodes;
  unsigned int g_spatial_total_exported_oocysts;
  unsigned int g_spatial_total_imported_human_infections;
  unsigned int g_spatial_total_imported_mosquito_infections;
  unsigned int g_spatial_imported_human_infection_counter;
  unsigned int g_spatial_imported_mosquito_infection_counter;
  unsigned int g_spatial_exported_barcode_counter;
  unsigned int g_spatial_exported_oocyst_counter;
  std::vector<boost::dynamic_bitset<> > g_spatial_exported_barcodes;
  std::vector<boost::dynamic_bitset<> > g_spatial_imported_barcodes;
  std::vector<boost::dynamic_bitset<> > g_spatial_exported_oocysts;
  std::vector<boost::dynamic_bitset<> > g_spatial_imported_oocysts;
  
  // island related
  static bool g_island_imports_plaf_linked_flag;
  unsigned int g_total_human_infections;
  unsigned int g_total_mosquito_infections;
  double g_percentage_imported_human_infections;
  double g_percentage_imported_mosquito_infections;
  
  // cotransmission and oocyst vectors
  std::vector<int> g_cotransmission_frequencies;
  unsigned int g_cotransmission_frequencies_counter;
  unsigned int g_cotransmission_frequencies_size;
  std::vector<int> g_oocyst_frequencies;
  unsigned int g_oocyst_frequencies_counter;
  unsigned int g_oocyst_frequencies_size;
  
  // demographic parameters;
  unsigned int g_N;
  int g_max_age;
  int g_average_age;
  // epidemiological parameters;
  double g_EIR;
  double g_a0;
  double g_rho;
  double g_zeta_meanlog;
  double g_zeta_sdlog;
  // diagnostic parameters;
  double g_ft;
  // entomological parameters;
  double g_mu0;
  double g_mean_mosquito_age;
  double g_beta_gradient;
  double g_beta_intercept;
  double g_ak;
  double g_Q0;
  // entomological time keeping variables
  int g_mosquito_deficit;
  unsigned int g_scourge_today;
  int g_mean_mv;
  std::vector<int> g_mosquito_next_biting_day_vector;
  int g_mosquito_biting_counter;
  int g_max_mosquito_biting_counter;
  std::vector<int> g_mosquito_eip_vector;
  int g_mosquito_eip_counter;
  int g_mosquito_current_eip;
  // delays and durations;
  double g_delay_mos;
  double g_delay_gam;
  double g_dur_E;
  double g_dur_T;
  double g_dur_D;
  double g_dur_U;
  double g_dur_P;
  double g_dur_A;
  double g_dur_AU;
  // Pre-erythrocytic immunity parameters;
  double g_d1;
  double g_dID;
  double g_ID0;
  double g_kD;
  double g_uD;
  double g_aD;
  double g_fD0;
  double g_gD;
  double g_alphaU;
  // Blood stage immunity parameters;
  double g_b0;
  double g_b1;
  double g_dB;
  double g_IB0;
  double g_kB;
  double g_uB;
  // Acquired and maternal immunity parameters;
  double g_phi0;
  double g_phi1;
  double g_dCA;
  double g_IC0;
  double g_kC;
  double g_uCA;
  double g_PM;
  double g_dCM;
  double g_mean_maternal_immunity;
  double g_sum_maternal_immunity;
  int g_total_mums;
  // Severe disease immunity parameters;
  double g_theta0;
  double g_theta1;
  double g_dVA;
  double g_IV0;
  double g_kV;
  double g_uVA;
  double g_fV0;
  double g_aV;
  double g_gammaV;
  double g_PVM;
  double g_dVM;
  double g_mean_maternal_severe_immunity;
  double g_sum_maternal_severe_immunity;
  // Newborns take maternal immunity from a random mother (aged 20 to <21 years) in their biting heterogeneity
  // group. g_het_brackets are the upper bounds of the groups (the equilibrium's het brackets). Mothers'
  // ICA and IVA are collected by group each day, and newborns draw from yesterday's mothers
  std::vector<double> g_het_brackets;
  std::vector<std::vector<double> > g_mother_ICA_today;
  std::vector<std::vector<double> > g_mother_IVA_today;
  std::vector<std::vector<double> > g_mother_ICA;
  std::vector<std::vector<double> > g_mother_IVA;
  // contributions to infectious reservoir by state and age;
  double g_gamma1;
  double g_cD;
  double g_cT;
  double g_cU;
  
  // barcode parameters
  static unsigned int g_identity_id;
  static unsigned int g_num_loci;
  static unsigned int g_ibd_length;
  static unsigned int g_barcode_length;
  static std::vector<double> g_plaf;
  static std::vector<double> g_prob_crossover;
  enum g_barcode_type_enum{
    ORDINARY, //0
    IBD, //1
    NUMBER_OF_BARCODE_TYPE_OPTIONS = 2
  };
  static g_barcode_type_enum g_barcode_type;
  
  // barcode drug related parameters
  bool g_resistance_flag;
  bool g_absolute_fitness_cost_flag;
  unsigned int g_number_of_resistance_loci;
  std::vector<unsigned int> g_resistance_loci;
  std::vector<double> g_cost_of_resistance;
  std::vector<unsigned int>  g_artemisinin_loci;
  std::vector<Drug> g_drugs;

  
  // mutation parameters
  bool g_mutation_flag;
  std::vector<double> g_mutation_rate;
  double g_mutation_treated_modifier;
  std::vector<unsigned int> g_mutations_today;
  unsigned int g_mutation_pos_allocator;
  
  // drug related parameters
  bool g_mft_flag;
  int g_drug_choice;
  unsigned int g_number_of_drugs;
  std::vector<double> g_partner_drug_ratios;
  double g_dur_SPC;
  
  // mosquito strain interaction params
  bool g_vector_adaptation_flag;
  std::vector<unsigned> g_vector_adaptation_loci;
  double g_local_oocyst_advantage;
  bool g_gametocyte_sterilisation_flag;
  double g_gametocyte_sterilisation;
  double g_oocyst_reduction_by_artemisinin;
  
  // non malaria fever parameters
  bool g_nmf_flag; // are we doing nmf work
  std::vector<double> g_mean_nmf_frequency; // mean number of days between nmf events
  std::vector<double> g_nmf_age_brackets; // probability that a non malarial fever happens
  double g_prob_of_testing_nmf;
  
  // Housekeeping parameters
  bool g_h_quiet_print;
  bool g_h_quiet_test_print;
  
  // Default Constructor
  Parameters();
  
  // Maternal immunity helpers
  static bool is_mother_age(int age) { return(age >= 20 * 365 && age < 21 * 365); }
  int het_group(double zeta) const;
  void set_default_het_brackets();
  void clear_mothers_today();
  void add_mother(double zeta, double ICA, double IVA);
  void start_mothers_day();
  bool draw_maternal_immunity(double zeta, double &ICA, double &IVA) const;
  
};


#endif