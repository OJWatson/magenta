//
//  magenta
//  main_saved_init.cpp
//
//  Created: OJ Watson on 06/12/2015
//
//  NOT ACCURATE - DO NOT USE CURRENTLY
//
//  Distributed under the MIT software licence
//
//  Further details (if any) of this set of functions can be found in the corresponding header file.
//
// ---------------------------------------------------------------------------

//#include <RcppArmadillo.h>
#include <iostream>
#include "parameters.h"
#include "probability.h"
#include <cassert> // for error checking
#include "person.h"
#include <chrono>
#include <functional>
#include <numeric>  
#include <algorithm>

using namespace std;
using namespace Rcpp;

// Create universe structure for all important variables
struct Universe {
  // Human storage
  std::vector<Person> population;
  std::vector<double> psi_vector;
  std::vector<double> zeta_vector;
  std::vector<double> pi_vector;
  // Mosquito storage
  std::vector<Mosquito> scourge;
  // Parameter storage
  Parameters parameters;
};


// ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// START: MAIN
// ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
//' Creates initial model simulation using a saved model state
//'
//' @param param_list parameter list generated with \code{Param_List_Simulation_Get_Create}
//' @return list with ptr to model state and loggers describing the current model state
//' @export
// [[Rcpp::export]]
Rcpp::List Simulation_Saved_Init_cpp(Rcpp::List param_list)
{
  
  // Initialise parameters
  Parameters parameters;

  auto int_matrix_to_bitsets = [](const std::vector<std::vector<int> > &matrix) {
    std::vector<boost::dynamic_bitset<> > res;
    res.reserve(matrix.size());
    for (const auto &row : matrix) {
      boost::dynamic_bitset<> bitset(Parameters::g_barcode_length);
      unsigned int limit = std::min<unsigned int>(row.size(), Parameters::g_barcode_length);
      for (unsigned int j = 0; j < limit; ++j) {
        bitset[j] = static_cast<bool>(row[j]);
      }
      res.emplace_back(bitset);
    }
    return res;
  };

  // Unpack R List to Rcpp Lists
  Rcpp::List savedState = param_list["savedState"];
  Rcpp::List population_List = savedState["population_List"];
  Rcpp::List populations_event_and_strains_List = savedState["populations_event_and_strains_List"];
  Rcpp::List scourge_List = savedState["scourge_List"];
  Rcpp::List parameters_List = savedState["parameters_List"];
  Rcpp::List parameter_housekeeping_List = savedState["parameter_housekeeping_List"];
  
  
  // prove that C++ code is being run
  parameters.g_h_quiet_print = Rcpp::as<bool>(parameter_housekeeping_List["g_h_quiet_print"]);
  parameters.g_h_quiet_test_print = Rcpp::as<bool>(parameter_housekeeping_List["g_h_quiet_test_print"]);
  rcpp_out(parameters.g_h_quiet_print, "Rcpp function is working!\n");
  
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // START: R -> C++ CONVERSIONS: parameters
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  // Update this first as needed for initialising the population
  parameters.g_N = Rcpp::as<unsigned int>(parameters_List["g_N"]);
  parameters.g_current_time = Rcpp::as<int>(parameters_List["g_current_time"]);
  parameters.g_mean_maternal_immunity = Rcpp::as<double>(parameters_List["g_mean_maternal_immunity"]);
  parameters.g_sum_maternal_immunity = Rcpp::as<double>(parameters_List["g_sum_maternal_immunity"]);
  parameters.g_total_mums = Rcpp::as<int>(parameters_List["g_total_mums"]);
  parameters.g_theta = Rcpp::as<vector<double> >(parameters_List["g_theta"]);
  parameters.g_calendar_day = Rcpp::as<int>(parameters_List["g_calendar_day"]);
  parameters.g_mosquito_deficit = Rcpp::as<int>(parameters_List["g_mosquito_deficit"]);
  parameters.g_scourge_today = Rcpp::as<int>(parameters_List["g_scourge_today"]);
  parameters.g_mean_mv = Rcpp::as<int>(parameters_List["g_mean_mv"]);

  if (parameters_List.containsElementNamed("g_years")) {
    parameters.g_years = Rcpp::as<double>(parameters_List["g_years"]);
  }
  if (parameters_List.containsElementNamed("g_max_age")) {
    parameters.g_max_age = Rcpp::as<int>(parameters_List["g_max_age"]);
  }
  if (parameters_List.containsElementNamed("g_average_age")) {
    parameters.g_average_age = Rcpp::as<int>(parameters_List["g_average_age"]);
  }
  if (parameters_List.containsElementNamed("g_EIR")) {
    parameters.g_EIR = Rcpp::as<double>(parameters_List["g_EIR"]);
  }
  if (parameters_List.containsElementNamed("g_a0")) {
    parameters.g_a0 = Rcpp::as<double>(parameters_List["g_a0"]);
  }
  if (parameters_List.containsElementNamed("g_rho")) {
    parameters.g_rho = Rcpp::as<double>(parameters_List["g_rho"]);
  }
  if (parameters_List.containsElementNamed("g_zeta_meanlog")) {
    parameters.g_zeta_meanlog = Rcpp::as<double>(parameters_List["g_zeta_meanlog"]);
  }
  if (parameters_List.containsElementNamed("g_zeta_sdlog")) {
    parameters.g_zeta_sdlog = Rcpp::as<double>(parameters_List["g_zeta_sdlog"]);
  }
  if (parameters_List.containsElementNamed("g_ft")) {
    parameters.g_ft = Rcpp::as<double>(parameters_List["g_ft"]);
  }
  if (parameters_List.containsElementNamed("g_mu0")) {
    parameters.g_mu0 = Rcpp::as<double>(parameters_List["g_mu0"]);
  }
  if (parameters_List.containsElementNamed("g_mean_mosquito_age")) {
    parameters.g_mean_mosquito_age = Rcpp::as<double>(parameters_List["g_mean_mosquito_age"]);
  }
  if (parameters_List.containsElementNamed("g_beta_gradient")) {
    parameters.g_beta_gradient = Rcpp::as<double>(parameters_List["g_beta_gradient"]);
  }
  if (parameters_List.containsElementNamed("g_beta_intercept")) {
    parameters.g_beta_intercept = Rcpp::as<double>(parameters_List["g_beta_intercept"]);
  }
  if (parameters_List.containsElementNamed("g_ak")) {
    parameters.g_ak = Rcpp::as<double>(parameters_List["g_ak"]);
  }
  if (parameters_List.containsElementNamed("g_Q0")) {
    parameters.g_Q0 = Rcpp::as<double>(parameters_List["g_Q0"]);
  }
  if (parameters_List.containsElementNamed("g_mosquito_next_biting_day_vector")) {
    parameters.g_mosquito_next_biting_day_vector = Rcpp::as<std::vector<int> >(parameters_List["g_mosquito_next_biting_day_vector"]);
  }
  if (parameters_List.containsElementNamed("g_mosquito_biting_counter")) {
    parameters.g_mosquito_biting_counter = Rcpp::as<int>(parameters_List["g_mosquito_biting_counter"]);
  }
  if (parameters_List.containsElementNamed("g_max_mosquito_biting_counter")) {
    parameters.g_max_mosquito_biting_counter = Rcpp::as<int>(parameters_List["g_max_mosquito_biting_counter"]);
  }
  if (parameters_List.containsElementNamed("g_delay_mos")) {
    parameters.g_delay_mos = Rcpp::as<double>(parameters_List["g_delay_mos"]);
  }
  if (parameters_List.containsElementNamed("g_delay_gam")) {
    parameters.g_delay_gam = Rcpp::as<double>(parameters_List["g_delay_gam"]);
  }
  if (parameters_List.containsElementNamed("g_dur_E")) {
    parameters.g_dur_E = Rcpp::as<double>(parameters_List["g_dur_E"]);
  }
  if (parameters_List.containsElementNamed("g_dur_T")) {
    parameters.g_dur_T = Rcpp::as<double>(parameters_List["g_dur_T"]);
  }
  if (parameters_List.containsElementNamed("g_dur_D")) {
    parameters.g_dur_D = Rcpp::as<double>(parameters_List["g_dur_D"]);
  }
  if (parameters_List.containsElementNamed("g_dur_U")) {
    parameters.g_dur_U = Rcpp::as<double>(parameters_List["g_dur_U"]);
  }
  if (parameters_List.containsElementNamed("g_dur_P")) {
    parameters.g_dur_P = Rcpp::as<double>(parameters_List["g_dur_P"]);
  }
  if (parameters_List.containsElementNamed("g_dur_A")) {
    parameters.g_dur_A = Rcpp::as<double>(parameters_List["g_dur_A"]);
  }
  if (parameters_List.containsElementNamed("g_dur_AU")) {
    parameters.g_dur_AU = Rcpp::as<double>(parameters_List["g_dur_AU"]);
  }
  if (parameters_List.containsElementNamed("g_d1")) {
    parameters.g_d1 = Rcpp::as<double>(parameters_List["g_d1"]);
  }
  if (parameters_List.containsElementNamed("g_dID")) {
    parameters.g_dID = Rcpp::as<double>(parameters_List["g_dID"]);
  }
  if (parameters_List.containsElementNamed("g_ID0")) {
    parameters.g_ID0 = Rcpp::as<double>(parameters_List["g_ID0"]);
  }
  if (parameters_List.containsElementNamed("g_kD")) {
    parameters.g_kD = Rcpp::as<double>(parameters_List["g_kD"]);
  }
  if (parameters_List.containsElementNamed("g_uD")) {
    parameters.g_uD = Rcpp::as<double>(parameters_List["g_uD"]);
  }
  if (parameters_List.containsElementNamed("g_aD")) {
    parameters.g_aD = Rcpp::as<double>(parameters_List["g_aD"]);
  }
  if (parameters_List.containsElementNamed("g_fD0")) {
    parameters.g_fD0 = Rcpp::as<double>(parameters_List["g_fD0"]);
  }
  if (parameters_List.containsElementNamed("g_gD")) {
    parameters.g_gD = Rcpp::as<double>(parameters_List["g_gD"]);
  }
  if (parameters_List.containsElementNamed("g_alphaU")) {
    parameters.g_alphaU = Rcpp::as<double>(parameters_List["g_alphaU"]);
  }
  if (parameters_List.containsElementNamed("g_b0")) {
    parameters.g_b0 = Rcpp::as<double>(parameters_List["g_b0"]);
  }
  if (parameters_List.containsElementNamed("g_b1")) {
    parameters.g_b1 = Rcpp::as<double>(parameters_List["g_b1"]);
  }
  if (parameters_List.containsElementNamed("g_dB")) {
    parameters.g_dB = Rcpp::as<double>(parameters_List["g_dB"]);
  }
  if (parameters_List.containsElementNamed("g_IB0")) {
    parameters.g_IB0 = Rcpp::as<double>(parameters_List["g_IB0"]);
  }
  if (parameters_List.containsElementNamed("g_kB")) {
    parameters.g_kB = Rcpp::as<double>(parameters_List["g_kB"]);
  }
  if (parameters_List.containsElementNamed("g_uB")) {
    parameters.g_uB = Rcpp::as<double>(parameters_List["g_uB"]);
  }
  if (parameters_List.containsElementNamed("g_phi0")) {
    parameters.g_phi0 = Rcpp::as<double>(parameters_List["g_phi0"]);
  }
  if (parameters_List.containsElementNamed("g_phi1")) {
    parameters.g_phi1 = Rcpp::as<double>(parameters_List["g_phi1"]);
  }
  if (parameters_List.containsElementNamed("g_dCA")) {
    parameters.g_dCA = Rcpp::as<double>(parameters_List["g_dCA"]);
  }
  if (parameters_List.containsElementNamed("g_IC0")) {
    parameters.g_IC0 = Rcpp::as<double>(parameters_List["g_IC0"]);
  }
  if (parameters_List.containsElementNamed("g_kC")) {
    parameters.g_kC = Rcpp::as<double>(parameters_List["g_kC"]);
  }
  if (parameters_List.containsElementNamed("g_uCA")) {
    parameters.g_uCA = Rcpp::as<double>(parameters_List["g_uCA"]);
  }
  if (parameters_List.containsElementNamed("g_PM")) {
    parameters.g_PM = Rcpp::as<double>(parameters_List["g_PM"]);
  }
  if (parameters_List.containsElementNamed("g_dCM")) {
    parameters.g_dCM = Rcpp::as<double>(parameters_List["g_dCM"]);
  }
  if (parameters_List.containsElementNamed("g_gamma1")) {
    parameters.g_gamma1 = Rcpp::as<double>(parameters_List["g_gamma1"]);
  }
  if (parameters_List.containsElementNamed("g_cD")) {
    parameters.g_cD = Rcpp::as<double>(parameters_List["g_cD"]);
  }
  if (parameters_List.containsElementNamed("g_cT")) {
    parameters.g_cT = Rcpp::as<double>(parameters_List["g_cT"]);
  }
  if (parameters_List.containsElementNamed("g_cU")) {
    parameters.g_cU = Rcpp::as<double>(parameters_List["g_cU"]);
  }
  if (parameters_List.containsElementNamed("g_total_human_infections")) {
    parameters.g_total_human_infections = Rcpp::as<unsigned int>(parameters_List["g_total_human_infections"]);
  }
  if (parameters_List.containsElementNamed("g_total_mosquito_infections")) {
    parameters.g_total_mosquito_infections = Rcpp::as<unsigned int>(parameters_List["g_total_mosquito_infections"]);
  }
  if (parameters_List.containsElementNamed("g_percentage_imported_human_infections")) {
    parameters.g_percentage_imported_human_infections = Rcpp::as<double>(parameters_List["g_percentage_imported_human_infections"]);
  }
  if (parameters_List.containsElementNamed("g_percentage_imported_mosquito_infections")) {
    parameters.g_percentage_imported_mosquito_infections = Rcpp::as<double>(parameters_List["g_percentage_imported_mosquito_infections"]);
  }
  if (parameters_List.containsElementNamed("g_cotransmission_frequencies")) {
    parameters.g_cotransmission_frequencies = Rcpp::as<std::vector<int> >(parameters_List["g_cotransmission_frequencies"]);
  }
  if (parameters_List.containsElementNamed("g_cotransmission_frequencies_counter")) {
    parameters.g_cotransmission_frequencies_counter = Rcpp::as<unsigned int>(parameters_List["g_cotransmission_frequencies_counter"]);
  }
  if (parameters_List.containsElementNamed("g_cotransmission_frequencies_size")) {
    parameters.g_cotransmission_frequencies_size = Rcpp::as<unsigned int>(parameters_List["g_cotransmission_frequencies_size"]);
  }
  if (parameters_List.containsElementNamed("g_oocyst_frequencies")) {
    parameters.g_oocyst_frequencies = Rcpp::as<std::vector<int> >(parameters_List["g_oocyst_frequencies"]);
  }
  if (parameters_List.containsElementNamed("g_oocyst_frequencies_counter")) {
    parameters.g_oocyst_frequencies_counter = Rcpp::as<unsigned int>(parameters_List["g_oocyst_frequencies_counter"]);
  }
  if (parameters_List.containsElementNamed("g_oocyst_frequencies_size")) {
    parameters.g_oocyst_frequencies_size = Rcpp::as<unsigned int>(parameters_List["g_oocyst_frequencies_size"]);
  }

  if (parameters_List.containsElementNamed("g_identity_id")) {
    Parameters::g_identity_id = Rcpp::as<unsigned int>(parameters_List["g_identity_id"]);
  }
  if (parameters_List.containsElementNamed("g_num_loci")) {
    Parameters::g_num_loci = Rcpp::as<unsigned int>(parameters_List["g_num_loci"]);
  }
  if (parameters_List.containsElementNamed("g_ibd_length")) {
    Parameters::g_ibd_length = Rcpp::as<unsigned int>(parameters_List["g_ibd_length"]);
  }
  if (parameters_List.containsElementNamed("g_barcode_length")) {
    Parameters::g_barcode_length = Rcpp::as<unsigned int>(parameters_List["g_barcode_length"]);
  }
  if (parameters_List.containsElementNamed("g_plaf")) {
    Parameters::g_plaf = Rcpp::as<std::vector<double> >(parameters_List["g_plaf"]);
  }
  if (parameters_List.containsElementNamed("g_prob_crossover")) {
    Parameters::g_prob_crossover = Rcpp::as<std::vector<double> >(parameters_List["g_prob_crossover"]);
  }
  if (parameters_List.containsElementNamed("g_barcode_type")) {
    Parameters::g_barcode_type = static_cast<Parameters::g_barcode_type_enum>(
      Rcpp::as<unsigned int>(parameters_List["g_barcode_type"])
    );
  }
  if (parameters_List.containsElementNamed("g_island_imports_plaf_linked_flag")) {
    Parameters::g_island_imports_plaf_linked_flag = Rcpp::as<bool>(parameters_List["g_island_imports_plaf_linked_flag"]);
  }

  parameters.g_spatial_type = static_cast<Parameters::g_spatial_type_enum>(
    Rcpp::as<unsigned int>(parameters_List["g_spatial_type"])
  );
  if (parameters_List.containsElementNamed("g_spatial_imported_cotransmission_frequencies")) {
    parameters.g_spatial_imported_cotransmission_frequencies = Rcpp::as<std::vector<int> >(parameters_List["g_spatial_imported_cotransmission_frequencies"]);
  }
  if (parameters_List.containsElementNamed("g_spatial_imported_oocyst_frequencies")) {
    parameters.g_spatial_imported_oocyst_frequencies = Rcpp::as<std::vector<int> >(parameters_List["g_spatial_imported_oocyst_frequencies"]);
  }
  if (parameters_List.containsElementNamed("g_spatial_exported_cotransmission_frequencies")) {
    parameters.g_spatial_exported_cotransmission_frequencies = Rcpp::as<std::vector<int> >(parameters_List["g_spatial_exported_cotransmission_frequencies"]);
  }
  if (parameters_List.containsElementNamed("g_spatial_exported_oocyst_frequencies")) {
    parameters.g_spatial_exported_oocyst_frequencies = Rcpp::as<std::vector<int> >(parameters_List["g_spatial_exported_oocyst_frequencies"]);
  }
  if (parameters_List.containsElementNamed("g_spatial_total_exported_barcodes")) {
    parameters.g_spatial_total_exported_barcodes = Rcpp::as<unsigned int>(parameters_List["g_spatial_total_exported_barcodes"]);
  }
  if (parameters_List.containsElementNamed("g_spatial_total_exported_oocysts")) {
    parameters.g_spatial_total_exported_oocysts = Rcpp::as<unsigned int>(parameters_List["g_spatial_total_exported_oocysts"]);
  }
  if (parameters_List.containsElementNamed("g_spatial_total_imported_human_infections")) {
    parameters.g_spatial_total_imported_human_infections = Rcpp::as<unsigned int>(parameters_List["g_spatial_total_imported_human_infections"]);
  }
  if (parameters_List.containsElementNamed("g_spatial_total_imported_mosquito_infections")) {
    parameters.g_spatial_total_imported_mosquito_infections = Rcpp::as<unsigned int>(parameters_List["g_spatial_total_imported_mosquito_infections"]);
  }
  if (parameters_List.containsElementNamed("g_spatial_imported_human_infection_counter")) {
    parameters.g_spatial_imported_human_infection_counter = Rcpp::as<unsigned int>(parameters_List["g_spatial_imported_human_infection_counter"]);
  }
  if (parameters_List.containsElementNamed("g_spatial_imported_mosquito_infection_counter")) {
    parameters.g_spatial_imported_mosquito_infection_counter = Rcpp::as<unsigned int>(parameters_List["g_spatial_imported_mosquito_infection_counter"]);
  }
  if (parameters_List.containsElementNamed("g_spatial_exported_barcode_counter")) {
    parameters.g_spatial_exported_barcode_counter = Rcpp::as<unsigned int>(parameters_List["g_spatial_exported_barcode_counter"]);
  }
  if (parameters_List.containsElementNamed("g_spatial_exported_oocyst_counter")) {
    parameters.g_spatial_exported_oocyst_counter = Rcpp::as<unsigned int>(parameters_List["g_spatial_exported_oocyst_counter"]);
  }
  if (parameters_List.containsElementNamed("g_spatial_exported_barcodes")) {
    parameters.g_spatial_exported_barcodes = int_matrix_to_bitsets(Rcpp::as<std::vector<std::vector<int> > >(parameters_List["g_spatial_exported_barcodes"]));
  }
  if (parameters_List.containsElementNamed("g_spatial_imported_barcodes")) {
    parameters.g_spatial_imported_barcodes = int_matrix_to_bitsets(Rcpp::as<std::vector<std::vector<int> > >(parameters_List["g_spatial_imported_barcodes"]));
  }
  if (parameters_List.containsElementNamed("g_spatial_exported_oocysts")) {
    parameters.g_spatial_exported_oocysts = int_matrix_to_bitsets(Rcpp::as<std::vector<std::vector<int> > >(parameters_List["g_spatial_exported_oocysts"]));
  }
  if (parameters_List.containsElementNamed("g_spatial_imported_oocysts")) {
    parameters.g_spatial_imported_oocysts = int_matrix_to_bitsets(Rcpp::as<std::vector<std::vector<int> > >(parameters_List["g_spatial_imported_oocysts"]));
  }

  if (parameters_List.containsElementNamed("g_resistance_flag")) {
    parameters.g_resistance_flag = Rcpp::as<bool>(parameters_List["g_resistance_flag"]);
  }
  if (parameters_List.containsElementNamed("g_absolute_fitness_cost_flag")) {
    parameters.g_absolute_fitness_cost_flag = Rcpp::as<bool>(parameters_List["g_absolute_fitness_cost_flag"]);
  }
  if (parameters_List.containsElementNamed("g_number_of_resistance_loci")) {
    parameters.g_number_of_resistance_loci = Rcpp::as<unsigned int>(parameters_List["g_number_of_resistance_loci"]);
  }
  if (parameters_List.containsElementNamed("g_resistance_loci")) {
    parameters.g_resistance_loci = Rcpp::as<std::vector<unsigned int> >(parameters_List["g_resistance_loci"]);
  }
  if (parameters_List.containsElementNamed("g_artemisinin_loci")) {
    parameters.g_artemisinin_loci = Rcpp::as<std::vector<unsigned int> >(parameters_List["g_artemisinin_loci"]);
  }
  if (parameters_List.containsElementNamed("g_cost_of_resistance")) {
    parameters.g_cost_of_resistance = Rcpp::as<std::vector<double> >(parameters_List["g_cost_of_resistance"]);
  }
  if (parameters_List.containsElementNamed("g_partner_drug_ratios")) {
    parameters.g_partner_drug_ratios = Rcpp::as<std::vector<double> >(parameters_List["g_partner_drug_ratios"]);
  }
  if (parameters_List.containsElementNamed("g_drug_choice")) {
    parameters.g_drug_choice = Rcpp::as<int>(parameters_List["g_drug_choice"]);
  }
  if (parameters_List.containsElementNamed("g_number_of_drugs")) {
    parameters.g_number_of_drugs = Rcpp::as<unsigned int>(parameters_List["g_number_of_drugs"]);
  }
  if (parameters_List.containsElementNamed("g_dur_SPC")) {
    parameters.g_dur_SPC = Rcpp::as<double>(parameters_List["g_dur_SPC"]);
  }
  parameters.g_drugs.clear();
  if (parameters_List.containsElementNamed("g_drugs")) {
    Rcpp::List drugs_list = parameters_List["g_drugs"];
    parameters.g_drugs.reserve(drugs_list.size());
    for (int i = 0; i < drugs_list.size(); ++i) {
      parameters.g_drugs.emplace_back(Rcpp::as<Rcpp::List>(drugs_list[i]));
    }
    parameters.g_number_of_drugs = static_cast<unsigned int>(parameters.g_drugs.size());
  }
  if (parameters_List.containsElementNamed("g_mft_flag")) {
    parameters.g_mft_flag = Rcpp::as<bool>(parameters_List["g_mft_flag"]);
  }

  if (parameters_List.containsElementNamed("g_mutation_flag")) {
    parameters.g_mutation_flag = Rcpp::as<bool>(parameters_List["g_mutation_flag"]);
  }
  if (parameters_List.containsElementNamed("g_mutation_rate")) {
    parameters.g_mutation_rate = Rcpp::as<std::vector<double> >(parameters_List["g_mutation_rate"]);
  }
  if (parameters_List.containsElementNamed("g_mutation_treated_modifier")) {
    parameters.g_mutation_treated_modifier = Rcpp::as<double>(parameters_List["g_mutation_treated_modifier"]);
  }
  if (parameters_List.containsElementNamed("g_mutations_today")) {
    parameters.g_mutations_today = Rcpp::as<std::vector<unsigned int> >(parameters_List["g_mutations_today"]);
  }
  if (parameters_List.containsElementNamed("g_mutation_pos_allocator")) {
    parameters.g_mutation_pos_allocator = Rcpp::as<unsigned int>(parameters_List["g_mutation_pos_allocator"]);
  }

  if (parameters_List.containsElementNamed("g_vector_adaptation_flag")) {
    parameters.g_vector_adaptation_flag = Rcpp::as<bool>(parameters_List["g_vector_adaptation_flag"]);
  }
  if (parameters_List.containsElementNamed("g_vector_adaptation_loci")) {
    parameters.g_vector_adaptation_loci = Rcpp::as<std::vector<unsigned int> >(parameters_List["g_vector_adaptation_loci"]);
  }
  if (parameters_List.containsElementNamed("g_local_oocyst_advantage")) {
    parameters.g_local_oocyst_advantage = Rcpp::as<double>(parameters_List["g_local_oocyst_advantage"]);
  }
  if (parameters_List.containsElementNamed("g_gametocyte_sterilisation_flag")) {
    parameters.g_gametocyte_sterilisation_flag = Rcpp::as<bool>(parameters_List["g_gametocyte_sterilisation_flag"]);
  }
  if (parameters_List.containsElementNamed("g_gametocyte_sterilisation")) {
    parameters.g_gametocyte_sterilisation = Rcpp::as<double>(parameters_List["g_gametocyte_sterilisation"]);
  }
  if (parameters_List.containsElementNamed("g_oocyst_reduction_by_artemisinin")) {
    parameters.g_oocyst_reduction_by_artemisinin = Rcpp::as<double>(parameters_List["g_oocyst_reduction_by_artemisinin"]);
  }

  if (parameters_List.containsElementNamed("g_nmf_flag")) {
    parameters.g_nmf_flag = Rcpp::as<bool>(parameters_List["g_nmf_flag"]);
  }
  if (parameters_List.containsElementNamed("g_mean_nmf_frequency")) {
    parameters.g_mean_nmf_frequency = Rcpp::as<std::vector<double> >(parameters_List["g_mean_nmf_frequency"]);
  }
  if (parameters_List.containsElementNamed("g_nmf_age_brackets")) {
    parameters.g_nmf_age_brackets = Rcpp::as<std::vector<double> >(parameters_List["g_nmf_age_brackets"]);
  }
  if (parameters_List.containsElementNamed("g_prob_of_testing_nmf")) {
    parameters.g_prob_of_testing_nmf = Rcpp::as<double>(parameters_List["g_prob_of_testing_nmf"]);
  }
  
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // END: R -> C++ CONVERSIONS: parameters
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  rcpp_out(parameters.g_h_quiet_print, "Paremeter list conversion is working!\n");
  
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // START: INITIALISATION
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  // Add the human population
  std::vector<Person> population;
  population.reserve(parameters.g_N);
  
  // Initialise vectors for keeping biting related variables
  std::vector<double> psi_vector(parameters.g_N);
  std::vector<double> zeta_vector(parameters.g_N);
  std::vector<double> pi_vector(parameters.g_N);
  
  // Add the mosquito population, i.e. scourge
  std::vector<Mosquito> scourge;
  unsigned int scourge_size = Rcpp::as<unsigned int>(scourge_List["Scourge_size"]);
  scourge.reserve(static_cast<int>(scourge_size));
  
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // START: R -> C++ CONVERSIONS: humans
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  std::vector<double> Zetas = Rcpp::as<vector<double> >(population_List["Zetas"]);
  std::vector<int> Infection_States_v = Rcpp::as<vector<int> >(population_List["Infection_States"]);
  std::vector<int> Ages_v = Rcpp::as<vector<int> >(population_List["Ages"]);
  std::vector<double> IB_v = Rcpp::as<vector<double> >(population_List["IB"]);
  std::vector<double> ICA_v = Rcpp::as<vector<double> >(population_List["ICA"]);
  std::vector<double> ICM_v = Rcpp::as<vector<double> >(population_List["ICM"]);
  std::vector<double> ID_v = Rcpp::as<vector<double> >(population_List["ID"]);
  std::vector<int> Treatment_Outcomes = Rcpp::as<vector<int> >(population_List["Treatment_Outcomes"]);
  std::vector<int> Recrudescence_Outcomes = Rcpp::as<vector<int> >(population_List["Recrudescence_Outcomes"]);
  std::vector<int> Drug_choices = Rcpp::as<vector<int> >(population_List["Drug_choices"]);
  std::vector<int> Slow_parasite_clearance = Rcpp::as<vector<int> >(population_List["Slow_parasite_clearance"]);
  std::vector<int> Day_of_nmf = Rcpp::as<vector<int> >(population_List["Day_of_nmf"]);
  std::vector<unsigned int> NMF_age_band = Rcpp::as<vector<unsigned int> >(population_List["NMF_age_band"]);
  std::vector<int> IB_last_boost_time = Rcpp::as<vector<int> >(population_List["IB_last_boost_time"]);
  std::vector<int> ICA_last_boost_time = Rcpp::as<vector<int> >(population_List["ICA_last_boost_time"]);
  std::vector<int> ID_last_boost_time = Rcpp::as<vector<int> >(population_List["ID_last_boost_time"]);
  std::vector<int> IB_last_calculated_time = Rcpp::as<vector<int> >(population_List["IB_last_calculated_time"]);
  std::vector<int> I_C_D_CM_last_calculated_time = Rcpp::as<vector<int> >(population_List["I_C_D_CM_last_calculated_time"]);
  std::vector<double> Immunity_boost_float = Rcpp::as<vector<double> >(population_List["Immunity_boost_float"]);
  
  std::vector<int> Day_of_InfectionStatus_change = Rcpp::as<vector<int> >(populations_event_and_strains_List["Day_of_InfectionStatus_change"]);
  std::vector<int> Day_of_strain_clearance = Rcpp::as<vector<int> >(populations_event_and_strains_List["Day_of_strain_clearance"]);
  std::vector<int> Day_of_death = Rcpp::as<vector<int> >(populations_event_and_strains_List["Day_of_death"]);
  std::vector<int> Day_of_last_treatment = Rcpp::as<vector<int> >(populations_event_and_strains_List["Day_of_last_treatment"]);
  std::vector<int> Day_of_prophylaxis_waning = Rcpp::as<vector<int> >(populations_event_and_strains_List["Day_of_prophylaxis_waning"]);
  std::vector<int> Number_of_Strains = Rcpp::as<vector<int> >(populations_event_and_strains_List["Number_of_Strains"]);
  std::vector<int> Number_of_Realised_Infections = Rcpp::as<vector<int> >(populations_event_and_strains_List["Number_of_Realised_Infections"]);
  
  Rcpp::List Infection_time_realisation_vectors = populations_event_and_strains_List["Infection_time_realisation_vectors"];
  Rcpp::List Infection_state_realisation_vectors = populations_event_and_strains_List["Infection_state_realisation_vectors"];
  Rcpp::List Infection_barcode_realisation_vectors = populations_event_and_strains_List["Infection_barcode_realisation_vectors"];
  
  Rcpp::List Strain_infection_state_vectors = populations_event_and_strains_List["Strain_infection_state_vectors"];
  Rcpp::List Strain_day_of_infection_state_change_vectors = populations_event_and_strains_List["Strain_day_of_infection_state_change_vectors"];
  Rcpp::List Strain_day_of_acquisition_vectors = populations_event_and_strains_List["Strain_day_of_acquisition_vectors"];
  Rcpp::List Strain_barcode_vectors = populations_event_and_strains_List["Strain_barcode_vectors"];
  
  
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // END: R -> C++ CONVERSIONS: humans
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  rcpp_out(parameters.g_h_quiet_print, "Human list conversions are working!\n");
  
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // START: R -> C++ CONVERSIONS: mosquitos
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  
  std::vector<int> Mosquito_Infection_States = Rcpp::as<vector<int> >(scourge_List["Mosquito_Infection_States"]);
  std::vector<bool> Mosquito_Off_Season = Rcpp::as<vector<bool> >(scourge_List["Mosquito_Off_Season"]);
  std::vector<unsigned int> Mosquito_Day_of_next_blood_meal = Rcpp::as<vector<unsigned int> >(scourge_List["Mosquito_Day_of_next_blood_meal"]);
  std::vector<unsigned int> Mosquito_Day_of_death = Rcpp::as<vector<unsigned int> >(scourge_List["Mosquito_Day_of_death"]);
  std::vector<unsigned int> Mosquito_Number_of_ruptured_oocysts = Rcpp::as<vector<unsigned int> >(scourge_List["Mosquito_Number_of_ruptured_oocysts"]);
  
  Rcpp::List Mosquito_Oocyst_rupture_time_vectors = scourge_List["Mosquito_Oocyst_rupture_time_vectors"];
  Rcpp::List Mosquito_Oocyst_remaining_spz_vectors = scourge_List["Mosquito_Oocyst_remaining_spz_vectors"];
  Rcpp::List Mosquito_Oocyst_barcode_male_vectors = scourge_List["Mosquito_Oocyst_barcode_male_vectors"];
  Rcpp::List Mosquito_Oocyst_barcode_female_vectors = scourge_List["Mosquito_Oocyst_barcode_female_vectors"];

  
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // END: R -> C++ CONVERSIONS: mosquitos
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  
  rcpp_out(parameters.g_h_quiet_print, "Mosquito list unpacking working!\n");
  
  // Use read in R equilibrium state to then allocate each individual accordingly given tehir age and biting heterogeneity
  
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // START: HUMAN FETCHING
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  // Necessary temporary veriables for unpacking strains within humans
  Strain temp_strain;
  std::vector<std::vector<bool> > temp_strain_barcode_vector;
  std::vector<int>  temp_strain_state_vector;
  std::vector<int>  temp_strain_state_change_time_vector;
  std::vector<int>  temp_strain_day_of_acquisition;
  unsigned int temp_barcode_iterator = 0;
  
  for (unsigned int n=0; n < parameters.g_N; n++) 
  {
    
    // Set their id
    population.emplace_back(parameters);
    population[n].set_m_person_ID(n);
    
    // Set infection state
    population[n].set_m_infection_state(static_cast<Person::InfectionStatus>(Infection_States_v[n]));
    
    // Set Age
    population[n].set_m_person_age(Ages_v[n]);
    
    // Set Immunities
    population[n].set_m_IB(IB_v[n]);
    population[n].set_m_ICA(ICA_v[n]);
    population[n].set_m_ICM(ICM_v[n]);
    population[n].set_m_ID(ID_v[n]);

    population[n].set_m_treatment_outcome(static_cast<Person::TreatmentOutcome>(Treatment_Outcomes[n]));
    population[n].set_m_recrudescence_outcome(static_cast<Person::RecrudescenceOutcome>(Recrudescence_Outcomes[n]));
    population[n].set_m_drug_choice(Drug_choices[n]);
    population[n].set_m_slow_parasite_clearance_bool(static_cast<bool>(Slow_parasite_clearance[n]));

    // Set Boost times
    population[n].set_m_IB_last_boost_time(IB_last_boost_time[n]);
    population[n].set_m_ICA_last_boost_time(ICA_last_boost_time[n]);
    population[n].set_m_ID_last_boost_time(ID_last_boost_time[n]);
    
    // Set Calc times
    population[n].set_m_IB_last_calculated_time(IB_last_calculated_time[n]);
    population[n].set_m_I_C_D_CM_last_calculated_time(I_C_D_CM_last_calculated_time[n]);
    
    // Set immunity boost float
    population[n].set_m_immunity_boost_float(Immunity_boost_float[n]);
    
    // Set Day Changes
    population[n].set_m_day_of_InfectionStatus_change(Day_of_InfectionStatus_change[n]);
    population[n].set_m_day_of_strain_clearance(Day_of_strain_clearance[n]);
    population[n].set_m_day_of_death(Day_of_death[n]);
    population[n].set_m_day_last_treated(Day_of_last_treatment[n]);
    population[n].set_m_day_prophylaxis_wanes(Day_of_prophylaxis_waning[n]);
    population[n].set_m_day_of_nmf_from_saved(Day_of_nmf[n]);
    population[n].set_m_nmf_age_band(NMF_age_band[n]);
    
    // Strain Numbers
    population[n].set_m_number_of_strains(Number_of_Strains[n]);
    
    // Strain Ralised Numbers
    population[n].set_m_number_of_realised_infections(Number_of_Realised_Infections[n]);
    
    // Set vectors
    population[n].set_m_infection_time_realisation_vector_from_vector(Rcpp::as<std::vector<int> >(Infection_time_realisation_vectors[n]));
    population[n].set_m_infection_state_realisation_vector_from_vector(Rcpp::as<std::vector<int> >(Infection_state_realisation_vectors[n]));
    population[n].set_m_infection_barcode_realisation_vector_from_vector_of_vector_bool(Rcpp::as<std::vector<std::vector<bool> > >(Infection_barcode_realisation_vectors[n]));
    
    // Set strains
    // -----------------------------------------
    
    // First grab the vector of times, states and barcode vectors for the indivdual
    temp_strain_state_vector = Rcpp::as<std::vector<int> >(Strain_infection_state_vectors[n]);
    temp_strain_state_change_time_vector = Rcpp::as<std::vector<int> >(Strain_day_of_infection_state_change_vectors[n]);
    temp_strain_day_of_acquisition = Rcpp::as<std::vector<int> >(Strain_day_of_acquisition_vectors[n]);
    temp_strain_barcode_vector = Rcpp::as<std::vector<std::vector<bool> > >(Strain_barcode_vectors[n]);
    
    // Loop over each vector and making a temporary strain which is to be pushed onto the human
    for(int s = 0; s < Number_of_Strains[n]; s++)
    {
    // Set the temp strains infection state and day of strain state change
    temp_strain.set_m_strain_infection_status(static_cast<Strain::InfectionStatus>(temp_strain_state_vector[s]));
    temp_strain.set_m_day_of_strain_infection_status_change(temp_strain_state_change_time_vector[s]);
    temp_strain.set_m_day_of_strain_acquisition(temp_strain_day_of_acquisition[s]);
    
    // fetch vector<bool> and turn into barcode and then add to strain
    for(temp_barcode_iterator = 0; temp_barcode_iterator < Parameters::g_barcode_length ; temp_barcode_iterator++ )
    {
      Strain::temp_barcode[temp_barcode_iterator] = temp_strain_barcode_vector[s][temp_barcode_iterator];
    }
    temp_strain.set_m_barcode(Strain::temp_barcode);
    
    // Push the created temp strain onto the person
    population[n].allocate_strain_with_push(temp_strain);
    }
    
    // -----------------------------------------
    
    // Set the next strain change day and next event day
    population[n].set_m_day_of_next_strain_state_change();
    population[n].set_m_day_of_next_event();
    
    // Asign zeta and add to the zeta vector which is required for calculating the overall probability of being bitten, pi
    population[n].set_m_individual_biting_rate(Zetas[n]);
    zeta_vector[n] = population[n].get_m_individual_biting_rate();
    
  }
  
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // END: HUMAN FETCHING
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  rcpp_out(parameters.g_h_quiet_print, "Human fetching working\n");
  
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // START: MOSQUITO FETCHING
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  rcpp_out(parameters.g_h_quiet_print, "Pre mosquito-fetching working!\n");
  int test = 0;
  for (unsigned int n=0; n < scourge_size; n++) 
  {
    
    // Set their id
    scourge.emplace_back(parameters);
    scourge[n].set_m_mosquito_ID(n);
    
    // Set infection state
    scourge[n].set_m_mosquito_infection_state(static_cast<Mosquito::InfectionStatus>(Mosquito_Infection_States[n]));
    if(scourge[n].get_m_mosquito_infection_state() == Mosquito::INFECTED){
      scourge[n].m_mosquito_infected = true;
    }
    
    // Set number of ruptured oocysts
    scourge[n].set_m_ruptured_oocyst_count(Mosquito_Number_of_ruptured_oocysts[n]);
    
    // Set Day Changes
    scourge[n].set_m_day_of_next_blood_meal(Mosquito_Day_of_next_blood_meal[n]);
    scourge[n].set_m_day_of_death(Mosquito_Day_of_death[n]);
    scourge[n].set_m_mosquito_off_season(Mosquito_Off_Season[n]);
    
    if(test ==0){
      rcpp_out(parameters.g_h_quiet_print, "Pre-mosquito_vectors working\n");
      test=1;
    }
    
    // Set vectors;
    scourge[n].set_m_oocyst_rupture_time_vector(Rcpp::as<std::vector<int> >(Mosquito_Oocyst_rupture_time_vectors[n]));
    scourge[n].set_m_oocyst_remaining_spz_count(Rcpp::as<std::vector<int> >(Mosquito_Oocyst_remaining_spz_vectors[n]));
    scourge[n].set_m_oocyst_barcode_male_vector_from_vector_of_vector_bool(Rcpp::as<std::vector<std::vector<bool> > >(Mosquito_Oocyst_barcode_male_vectors[n]));
    scourge[n].set_m_oocyst_barcode_female_vector_from_vector_of_vector_bool(Rcpp::as<std::vector<std::vector<bool> > >(Mosquito_Oocyst_barcode_female_vectors[n]));
    
    // Set the next event day
    scourge[n].schedule_m_day_of_next_event();

    
  }
  
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // END: MOSQUITO FETCHING
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  rcpp_out(parameters.g_h_quiet_print, "Mosquito fetching working\n");
  
  /*
  
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // END: INITIALISATION
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  */
  
  // Final infection states
  std::vector<int> Infection_States(parameters.g_N);
  std::vector<double> Ages(parameters.g_N);
  std::vector<double> IB(parameters.g_N);
  std::vector<double> ICA(parameters.g_N);
  std::vector<double> ICM(parameters.g_N);
  std::vector<double> ID(parameters.g_N);
  
  // status eq for logging and other logging variables
  std::vector<double> status_eq = { 0,0,0,0,0,0 };
  
  // loop through population and grab ages and immunities for error checking
  for (unsigned int element = 0; element < parameters.g_N ; element++) 
  {
    
    // Match infection state and schedule associated next state change
    switch (population[element].get_m_infection_state())
    {
    case Person::SUSCEPTIBLE:
      status_eq[0]++;
      break;
    case Person::DISEASED:
      status_eq[1]++;
      break;
    case Person::ASYMPTOMATIC:
      status_eq[2]++;
      break;
    case Person::SUBPATENT:
      status_eq[3]++;
      break;
    case Person::TREATED:
      status_eq[4]++;
      break;
    case Person::PROPHYLAXIS:
      status_eq[5]++;
      break;
    default:
      assert(NULL && "Schedule Infection Status Change Error - person's infection status not S, D, A, U, T or P");
    break;
    }
    
    // Ages and immunity 
    // TODO: Figure out the best way of standardising this logging 
    // Something like passing in a function name within the param_list which is the 
    // name for a logger written else where which then returns the Loggers obeject below
    Infection_States[element] = static_cast<int>(population[element].get_m_infection_state());
    population[element].update_immunities_to_today(parameters);
    Ages[element] = population[element].get_m_person_age();
    IB[element] = population[element].get_m_IB();
    ICA[element] = population[element].get_m_ICA();
    ICM[element] = population[element].get_m_ICM();
    ID[element] = population[element].get_m_ID();
    
  }
  
  // divide by population size and log counter and print to give overview
  rcpp_out(parameters.g_h_quiet_print, "S | D | A | U | T | P:\n");
  
  for (int element = 0; element < 6; element++) 
  {
    status_eq[element] /= (parameters.g_N);
    rcpp_out(parameters.g_h_quiet_print, std::to_string(status_eq[element]) + " | ");
  }
  
  // Create Rcpp loggers list
  Rcpp::List Loggers = Rcpp::List::create(
    Rcpp::Named("S")=status_eq[0] ,
    Rcpp::Named("D")=status_eq[1],
    Rcpp::Named("A")=status_eq[2],
    Rcpp::Named("U")=status_eq[3],
    Rcpp::Named("T")=status_eq[4],
    Rcpp::Named("P")=status_eq[5],
    Rcpp::Named("InfectionStates")=Infection_States, 
    Rcpp::Named("Ages")=Ages, 
    Rcpp::Named("IB")=IB,
    Rcpp::Named("ICA")=ICA,
    Rcpp::Named("ICM")=ICM,
    Rcpp::Named("ID")=ID
  );
  
  
  // Create universe ptr for memory-continuiation
  Rcpp::XPtr<Universe> universe_ptr(new Universe{ population, psi_vector, zeta_vector, pi_vector, scourge, parameters},
                                    true);
  
  // Return Named List with pointer and loggers
  return Rcpp::List::create(Rcpp::Named("Ptr") = universe_ptr, Rcpp::Named("Loggers")=Loggers);
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // fini
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
}


