//
//  magenta
//  main_get.cpp
//
//  Created: OJ Watson on 06/12/2015
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
//#define NDEBUG; // This causes all asserts to not be used - good for use after debugging and development has finished

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
//' Returns whole model to R in series of nested lists
//'
//' @param param_list parameter list generated with \code{Param_List_Simulation_Get_Create}
//' @return list of 4 lists with the entire model state
//' @export
// [[Rcpp::export]]
Rcpp::List Simulation_Get_cpp(Rcpp::List param_list)
{
  

  // start timer
  chrono::high_resolution_clock::time_point t0 = std::chrono::high_resolution_clock::now();
  
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // START: R -> C++ CONVERSIONS
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  // Create universe pointer from param_list statePtr
  Rcpp::XPtr<Universe> universe_ptr = Rcpp::as<Rcpp::XPtr<Universe> > (param_list["statePtr"]);
  
    // prove that C++ code is being run
  rcpp_out(universe_ptr->parameters.g_h_quiet_print, "Rcpp function is working!\n");
  
  // Initialise variables from the statePtr provided that are needing to be saved
  std::vector<double> Zeta = universe_ptr->zeta_vector;
  
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // END: R -> C++ CONVERSIONS
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  rcpp_out(universe_ptr->parameters.g_h_quiet_print, "Pointer unpacking working!\n");
  
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // START: MODEL STATE GET
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // START: HUMAN POPULATION
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  // Create model state storage
  std::vector<int> Infection_States(universe_ptr->parameters.g_N);
  std::vector<int> Ages(universe_ptr->parameters.g_N);
  std::vector<double> IB(universe_ptr->parameters.g_N);
  std::vector<double> ICA(universe_ptr->parameters.g_N);
  std::vector<double> ICM(universe_ptr->parameters.g_N);
  std::vector<double> cA(universe_ptr->parameters.g_N);
  std::vector<double> ID(universe_ptr->parameters.g_N);
  std::vector<int> Treatment_Outcomes(universe_ptr->parameters.g_N);
  std::vector<int> Recrudescence_Outcomes(universe_ptr->parameters.g_N);
  std::vector<int> Drug_choices(universe_ptr->parameters.g_N);
  std::vector<int> Slow_parasite_clearance(universe_ptr->parameters.g_N);
  std::vector<int> Day_of_nmf(universe_ptr->parameters.g_N);
  std::vector<unsigned int> NMF_age_band(universe_ptr->parameters.g_N);
  std::vector<double> IB_last_boost_time(universe_ptr->parameters.g_N);
  std::vector<double> ICA_last_boost_time(universe_ptr->parameters.g_N);
  std::vector<double> ID_last_boost_time(universe_ptr->parameters.g_N);
  std::vector<int> IB_last_calculated_time(universe_ptr->parameters.g_N);
  std::vector<int> I_C_D_CM_last_calculated_time(universe_ptr->parameters.g_N);
  std::vector<double> Immunity_boost_float(universe_ptr->parameters.g_N);
  std::vector<int>  Day_of_InfectionStatus_change(universe_ptr->parameters.g_N);	
  std::vector<int>  Day_of_strain_clearance(universe_ptr->parameters.g_N);		
  std::vector<int>  Day_of_death(universe_ptr->parameters.g_N);					
  std::vector<int>  Number_of_Strains(universe_ptr->parameters.g_N);	
  std::vector<int> Number_of_Realised_Infections(universe_ptr->parameters.g_N);
  std::vector<int> Day_of_next_strain_state_change(universe_ptr->parameters.g_N);
  std::vector<int> Day_of_next_event(universe_ptr->parameters.g_N);
  std::vector<int> Day_of_last_treatment(universe_ptr->parameters.g_N);
  std::vector<int> Day_of_prophylaxis_waning(universe_ptr->parameters.g_N);
  std::vector<std::vector<int> > Infection_time_realisation_vectors(universe_ptr->parameters.g_N);
  std::vector<std::vector<int> > Infection_state_realisation_vectors(universe_ptr->parameters.g_N);
  std::vector<std::vector<std::vector<bool> > > Infection_barcode_realisation_vectors(universe_ptr->parameters.g_N);
  std::vector<std::vector<std::vector<bool> > > Strain_barcode_vectors(universe_ptr->parameters.g_N);
  std::vector<std::vector<int> > Strain_infection_state_vectors(universe_ptr->parameters.g_N);
  std::vector<std::vector<int> > Strain_day_of_infection_state_change_vectors(universe_ptr->parameters.g_N);
  std::vector<std::vector<int> > Strain_day_of_acquisition_vectors(universe_ptr->parameters.g_N);
  std::vector<std::vector<bool> > Strain_cotransmission(universe_ptr->parameters.g_N);
  
  std::vector<std::vector<unsigned int> > recent_barcode_integers(universe_ptr->parameters.g_N);
  
  rcpp_out(universe_ptr->parameters.g_h_quiet_print, "Vector initialisation working!\n");
  
  // Temporary necessities for casting vectors for pending states
  std::vector<Person::InfectionStatus> temp_infection_state_realisation_vector{};
  std::vector<bool> temp_cotransmission_vector{};
  unsigned int temp_status_iterator = 0;
  
  // Temporary necessities for pending barcodes
  std::vector<boost::dynamic_bitset<> > temp_infection_barcode_realisation_vector = {};
  std::vector<bool> temp_barcode_bool_vector(universe_ptr->parameters.g_barcode_length,false);
  unsigned int temp_barcode_iterator = 0;
  
  // Temporary necessities for strains
  Strain temp_strain;
  boost::dynamic_bitset<> temp_barcode;
  int temp_strain_iterator = 0;

  rcpp_out(universe_ptr->parameters.g_h_quiet_print, "Preloop working!\n");

  
  for (unsigned int element = 0; element < universe_ptr->parameters.g_N ; element++) 
  {
    
    // Infection States
    Infection_States[element] = static_cast<int>(universe_ptr->population[element].get_m_infection_state());
    
    // Ages
    Ages[element] = universe_ptr->population[element].get_m_person_age();
    
    // Immunities
    IB[element] = universe_ptr->population[element].get_m_IB();
    ICA[element] = universe_ptr->population[element].get_m_ICA();
    ICM[element] = universe_ptr->population[element].get_m_ICM();
    ID[element] = universe_ptr->population[element].get_m_ID();
    cA[element] = universe_ptr->population[element].get_m_cA();
    Treatment_Outcomes[element] = static_cast<int>(universe_ptr->population[element].get_m_treatment_outcome());
    Recrudescence_Outcomes[element] = static_cast<int>(universe_ptr->population[element].get_m_recrudescence_outcome());
    Drug_choices[element] = universe_ptr->population[element].get_m_drug_choice();
    Slow_parasite_clearance[element] = static_cast<int>(universe_ptr->population[element].get_m_slow_parasite_clearance_bool());
    Day_of_nmf[element] = universe_ptr->population[element].get_m_day_of_nmf();
    NMF_age_band[element] = universe_ptr->population[element].get_m_nmf_age_band();
    // Boost times
    IB_last_boost_time[element] = universe_ptr->population[element].get_m_IB_last_boost_time();
    ICA_last_boost_time[element] = universe_ptr->population[element].get_m_ICA_last_boost_time();
    ID_last_boost_time[element] = universe_ptr->population[element].get_m_ID_last_boost_time();
    
    // Calc times
    IB_last_calculated_time[element] = universe_ptr->population[element].get_m_IB_last_calculated_time();
    I_C_D_CM_last_calculated_time[element] = universe_ptr->population[element].get_m_I_C_D_CM_last_calculated_time();
    
    // Immunity float
    Immunity_boost_float[element] = universe_ptr->population[element].get_m_immunity_boost_float();
    
    // Day Changes
    Day_of_InfectionStatus_change[element] = universe_ptr->population[element].get_m_day_of_InfectionStatus_change();
    Day_of_strain_clearance[element] = universe_ptr->population[element].get_m_day_of_strain_clearance();
    Day_of_death[element]	 = universe_ptr->population[element].get_m_day_of_death();
    Day_of_last_treatment[element] = universe_ptr->population[element].get_m_day_last_treated(); 
    Day_of_prophylaxis_waning[element] = universe_ptr->population[element].get_m_day_prophylaxis_wanes(); 
    
    // Strain Numbers
    Number_of_Strains[element] = universe_ptr->population[element].get_m_number_of_strains();
    
    // Temp strain to next change state
    Day_of_next_strain_state_change[element] = universe_ptr->population[element].get_m_day_of_next_strain_state_change();
    Day_of_next_event[element] = universe_ptr->population[element].get_m_day_of_next_event();
    
    // Realised Infections
    Number_of_Realised_Infections[element] = universe_ptr->population[element].get_m_number_of_realised_infections();
    
    // Pending Infection time vector
    Infection_time_realisation_vectors[element] = universe_ptr->population[element].get_m_infection_time_realisation_vector();
    
     rcpp_out(universe_ptr->parameters.g_h_quiet_print, "Prepending working!" + std::to_string(element) + "\n");
    
    // Pending Infection barcode and state vector
    // ---------------------------------------
    temp_infection_barcode_realisation_vector = universe_ptr->population[element].get_m_infection_barcode_realisation_vector();
    Infection_barcode_realisation_vectors[element].reserve(temp_infection_barcode_realisation_vector.size());
    
    temp_infection_state_realisation_vector = universe_ptr->population[element].get_m_infection_state_realisation_vector();
    Infection_state_realisation_vectors[element].reserve(temp_infection_state_realisation_vector.size());
    
    rcpp_out(universe_ptr->parameters.g_h_quiet_print, "Prestatus  working!" + std::to_string(element) + "\n");
    for(temp_status_iterator = 0 ; temp_status_iterator < static_cast<unsigned int>(temp_infection_state_realisation_vector.size()) ; )
    {
      Infection_state_realisation_vectors[element].emplace_back(temp_infection_state_realisation_vector[temp_status_iterator]);

      for(temp_barcode_iterator = 0; temp_barcode_iterator < universe_ptr->parameters.g_barcode_length ; temp_barcode_iterator++ )
      {
        temp_barcode_bool_vector[temp_barcode_iterator] = static_cast<bool>(temp_infection_barcode_realisation_vector[temp_status_iterator][temp_barcode_iterator]);
      }

      Infection_barcode_realisation_vectors[element].push_back(temp_barcode_bool_vector);
      temp_status_iterator++;
    }
    
    temp_infection_barcode_realisation_vector.clear();
    temp_infection_state_realisation_vector.clear();

    // Active strains
    // ---------------------------------------
    // Reserve space for active strains
    Strain_barcode_vectors[element].reserve(Number_of_Strains[element]);
    Strain_infection_state_vectors[element].reserve(Number_of_Strains[element]);
    Strain_day_of_infection_state_change_vectors[element].reserve(Number_of_Strains[element]);
    
    std::vector<Strain> temp_strain_vector;
    temp_strain_vector = universe_ptr->population[element].get_m_active_strains();
    
    // rcpp_out(universe_ptr->parameters.g_h_quiet_print, "Prestrains loop working!\n");
    
    // Loop through each strain converting into barcodes, states and days of infection state changes
    for(temp_strain_iterator = 0 ; temp_strain_iterator < Number_of_Strains[element] ; temp_strain_iterator++)
    {
      
      temp_strain = universe_ptr->population[element].get_m_person_strain_x(temp_strain_iterator);
      
      Strain_infection_state_vectors[element].push_back(temp_strain.get_m_strain_infection_status());
      Strain_day_of_infection_state_change_vectors[element].push_back(temp_strain.get_m_day_of_strain_infection_status_change());
      Strain_day_of_acquisition_vectors[element].push_back(temp_strain.get_m_day_of_strain_acquisition());
      // fetch barcode and turn into vector<bool>
      temp_barcode = temp_strain.get_m_barcode();
      for(temp_barcode_iterator = 0; temp_barcode_iterator < universe_ptr->parameters.g_barcode_length ; temp_barcode_iterator++ )
      {
        temp_barcode_bool_vector[temp_barcode_iterator] = temp_barcode[temp_barcode_iterator];
      }
      Strain_barcode_vectors[element].push_back(temp_barcode_bool_vector);
      
      if(temp_strain_iterator == (Number_of_Strains[element]-1)){
        if(universe_ptr->parameters.g_barcode_type == Parameters::IBD){
          recent_barcode_integers[element].reserve(1);
          recent_barcode_integers[element] = Strain::ibd_barcode_to_integer_vector(temp_barcode);
        }
      }
    }
    
     rcpp_out(universe_ptr->parameters.g_h_quiet_print, "Poststrains loop working!\n");
    
  }
    
    // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    // END: HUMAN POPULATION
    // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    
    // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    // START: MOSQUITO POPULATION
    // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    
    // Create model state storage
    unsigned int scourge_size = universe_ptr->scourge.size();
  
    std::vector<int> Mosquito_Infection_States(scourge_size);
    std::vector<int> Mosquito_Day_of_next_blood_meal(scourge_size);	
    std::vector<bool> Mosquito_Off_Season(scourge_size);	
    std::vector<int> Mosquito_Day_of_death(scourge_size);
    std::vector<int> Mosquito_Number_of_ruptured_oocysts(scourge_size);
    std::vector<std::vector<int> > Mosquito_Oocyst_remaining_spz_vectors(scourge_size);
    std::vector<std::vector<int> > Mosquito_Oocyst_rupture_time_vectors(scourge_size);
    std::vector<std::vector<std::vector<bool> > > Mosquito_Oocyst_barcode_male_vectors(scourge_size);
    std::vector<std::vector<std::vector<bool> > > Mosquito_Oocyst_barcode_female_vectors(scourge_size);
    
    // Temporary necessities for pending barcodes
    std::vector<boost::dynamic_bitset<>> temp_male_barcode_realisation_vector = {};
    std::vector<boost::dynamic_bitset<>> temp_female_barcode_realisation_vector = {};

    std::vector<bool> temp_barcode_male_bool_vector{};
    std::vector<bool> temp_barcode_female_bool_vector{};
    temp_barcode_female_bool_vector.reserve(universe_ptr->parameters.g_barcode_length);
    temp_barcode_male_bool_vector.reserve(universe_ptr->parameters.g_barcode_length);
    
    rcpp_out(universe_ptr->parameters.g_h_quiet_print, "Premosquito loop working!\n");

    for (unsigned int element = 0; element < scourge_size ; element++) 
    {
      
      // Infection States
      Mosquito_Infection_States[element] = static_cast<int>(universe_ptr->scourge[element].get_m_mosquito_infection_state());

      // Day Changes
      Mosquito_Day_of_death[element] = universe_ptr->scourge[element].get_m_day_of_death();

      Mosquito_Day_of_next_blood_meal[element] = universe_ptr->scourge[element].get_m_day_of_next_blood_meal();

      Mosquito_Off_Season[element] = universe_ptr->scourge[element].get_m_mosquito_off_season();
      
      // Ruptured oocyst Numbers 
      Mosquito_Number_of_ruptured_oocysts[element] = static_cast<int>(universe_ptr->scourge[element].get_m_ruptured_oocyst_count());

      // Pending Rupture time vector
      Mosquito_Oocyst_rupture_time_vectors[element] = universe_ptr->scourge[element].get_m_oocyst_rupture_time_vector();
      
      // Pending Infection time vector
      Mosquito_Oocyst_remaining_spz_vectors[element] = universe_ptr->scourge[element].get_m_oocyst_remaining_spz_count();

      // Pending oocyst barcode vector
      // ---------------------------------------
      temp_male_barcode_realisation_vector = universe_ptr->scourge[element].get_m_oocyst_barcode_male_vector();

      temp_female_barcode_realisation_vector = universe_ptr->scourge[element].get_m_oocyst_barcode_female_vector();

      Mosquito_Oocyst_barcode_male_vectors[element].reserve(temp_male_barcode_realisation_vector.size());

      Mosquito_Oocyst_barcode_female_vectors[element].reserve(temp_female_barcode_realisation_vector.size());

      temp_status_iterator = 0;

      for(temp_status_iterator = 0; temp_status_iterator < temp_male_barcode_realisation_vector.size() ; temp_status_iterator++)
      {

        // fetch barcode and turn into vector<bool>
        for(temp_barcode_iterator = 0; temp_barcode_iterator < universe_ptr->parameters.g_barcode_length ; temp_barcode_iterator++ )
        {
          temp_barcode_male_bool_vector.push_back(temp_male_barcode_realisation_vector[temp_status_iterator][temp_barcode_iterator]);
          temp_barcode_female_bool_vector.push_back(temp_female_barcode_realisation_vector[temp_status_iterator][temp_barcode_iterator]);
        }
        
        Mosquito_Oocyst_barcode_male_vectors[element].push_back(temp_barcode_male_bool_vector);
        Mosquito_Oocyst_barcode_female_vectors[element].push_back(temp_barcode_female_bool_vector);
        temp_barcode_male_bool_vector.clear();
        temp_barcode_female_bool_vector.clear();

      }

      temp_male_barcode_realisation_vector.clear();
      temp_female_barcode_realisation_vector.clear();

      // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
      // END: MOSQUITO POPULATION
      // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    }
    
  
  
  rcpp_out(universe_ptr->parameters.g_h_quiet_print, "Postloop working!\n");
    
  
  // Create Rcpp population list
  Rcpp::List population_List = Rcpp::List::create(
    Rcpp::Named("Infection_States")=Infection_States,
    Rcpp::Named("Zetas")=Zeta,
    Rcpp::Named("Ages")=Ages,
    Rcpp::Named("IB")=IB,
    Rcpp::Named("ICA")=ICA,
    Rcpp::Named("ICM")=ICM,
    Rcpp::Named("ID")=ID,
    Rcpp::Named("cA")=cA,
    Rcpp::Named("Treatment_Outcomes")=Treatment_Outcomes,
    Rcpp::Named("Recrudescence_Outcomes")=Recrudescence_Outcomes,
    Rcpp::Named("Drug_choices")=Drug_choices,
    Rcpp::Named("Slow_parasite_clearance")=Slow_parasite_clearance,
    Rcpp::Named("Day_of_nmf")=Day_of_nmf,
    Rcpp::Named("NMF_age_band")=NMF_age_band,
    Rcpp::Named("IB_last_boost_time")=IB_last_boost_time,
    Rcpp::Named("ICA_last_boost_time")=ICA_last_boost_time,
    Rcpp::Named("ID_last_boost_time")=ID_last_boost_time,
    Rcpp::Named("IB_last_calculated_time")=IB_last_calculated_time,
    Rcpp::Named("I_C_D_CM_last_calculated_time")=I_C_D_CM_last_calculated_time,
    Rcpp::Named("Immunity_boost_float")=Immunity_boost_float
  );
  
  // Create Rcpp population event_and_strain list
  Rcpp::List populations_event_and_strains_List = Rcpp::List::create(
    Rcpp::Named("Day_of_InfectionStatus_change")=Day_of_InfectionStatus_change,
    Rcpp::Named("Day_of_strain_clearance")=Day_of_strain_clearance,
    Rcpp::Named("Day_of_death")=Day_of_death,
    Rcpp::Named("Number_of_Strains")=Number_of_Strains,
    Rcpp::Named("Day_of_next_strain_state_change")=Day_of_next_strain_state_change,
    Rcpp::Named("Day_of_next_event")=Day_of_next_event,
    Rcpp::Named("Day_of_last_treatment")=Day_of_last_treatment,
    Rcpp::Named("Day_of_prophylaxis_waning")=Day_of_prophylaxis_waning,
    Rcpp::Named("Number_of_Realised_Infections")=Number_of_Realised_Infections,
    Rcpp::Named("Infection_time_realisation_vectors")=Infection_time_realisation_vectors,
    Rcpp::Named("Infection_state_realisation_vectors")=Infection_state_realisation_vectors,
    Rcpp::Named("Infection_barcode_realisation_vectors")=Infection_barcode_realisation_vectors,
    Rcpp::Named("Strain_infection_state_vectors")=Strain_infection_state_vectors,
    Rcpp::Named("Strain_day_of_infection_state_change_vectors")=Strain_day_of_infection_state_change_vectors,
    Rcpp::Named("Strain_day_of_acquisition_vectors")=Strain_day_of_acquisition_vectors,
    Rcpp::Named("Strain_barcode_vectors")=Strain_barcode_vectors,
    Rcpp::Named("Strain_cotransmission")=Strain_cotransmission,
    Rcpp::Named("Recent_identity_vectors")=recent_barcode_integers
  );
  
  // Create Mosquito population list
  Rcpp::List scourge_List = Rcpp::List::create(
    Rcpp::Named("Mosquito_Infection_States")=Mosquito_Infection_States,
    Rcpp::Named("Mosquito_Day_of_next_blood_meal")=Mosquito_Day_of_next_blood_meal,
    Rcpp::Named("Mosquito_Off_Season")=Mosquito_Off_Season,
    Rcpp::Named("Mosquito_Day_of_death")=Mosquito_Day_of_death,
    Rcpp::Named("Mosquito_Number_of_ruptured_oocysts")=Mosquito_Number_of_ruptured_oocysts,
    Rcpp::Named("Mosquito_Oocyst_remaining_spz_vectors")=Mosquito_Oocyst_remaining_spz_vectors,
    Rcpp::Named("Mosquito_Oocyst_rupture_time_vectors")=Mosquito_Oocyst_rupture_time_vectors,
    Rcpp::Named("Mosquito_Oocyst_barcode_male_vectors")=Mosquito_Oocyst_barcode_male_vectors,
    Rcpp::Named("Mosquito_Oocyst_barcode_female_vectors")=Mosquito_Oocyst_barcode_female_vectors,
    Rcpp::Named("Scourge_size")=universe_ptr->scourge.size()
  );
  
  // Create Rcpp Parameters list

  // drugs list
  std::vector<Rcpp::List> drugs;
  drugs.reserve(universe_ptr->parameters.g_number_of_drugs);
  for (auto d : universe_ptr->parameters.g_drugs) {
    drugs.emplace_back(d.drug_to_rcpp_list());
  }

  auto bitset_matrix = [&](const std::vector<boost::dynamic_bitset<> > &bitsets) {
    std::vector<std::vector<int> > res(bitsets.size(), std::vector<int>(Parameters::g_barcode_length));
    for (size_t i = 0; i < bitsets.size(); ++i) {
      for (unsigned int j = 0; j < Parameters::g_barcode_length; ++j) {
        res[i][j] = static_cast<int>(bitsets[i][j]);
      }
    }
    return res;
  };

  std::vector<std::vector<int> > spatial_exported_barcodes;
  if (!universe_ptr->parameters.g_spatial_exported_barcodes.empty()) {
    spatial_exported_barcodes = bitset_matrix(universe_ptr->parameters.g_spatial_exported_barcodes);
  }

  std::vector<std::vector<int> > spatial_imported_barcodes;
  if (!universe_ptr->parameters.g_spatial_imported_barcodes.empty()) {
    spatial_imported_barcodes = bitset_matrix(universe_ptr->parameters.g_spatial_imported_barcodes);
  }

  std::vector<std::vector<int> > spatial_exported_oocysts;
  if (!universe_ptr->parameters.g_spatial_exported_oocysts.empty()) {
    spatial_exported_oocysts = bitset_matrix(universe_ptr->parameters.g_spatial_exported_oocysts);
  }

  std::vector<std::vector<int> > spatial_imported_oocysts;
  if (!universe_ptr->parameters.g_spatial_imported_oocysts.empty()) {
    spatial_imported_oocysts = bitset_matrix(universe_ptr->parameters.g_spatial_imported_oocysts);
  }

  Rcpp::List parameters_List = Rcpp::List::create(
    Rcpp::Named("g_current_time")=universe_ptr->parameters.g_current_time,
    Rcpp::Named("g_calendar_day")=universe_ptr->parameters.g_calendar_day,
    Rcpp::Named("g_theta")=universe_ptr->parameters.g_theta,
    Rcpp::Named("g_years")=universe_ptr->parameters.g_years,
    Rcpp::Named("g_mean_maternal_immunity")=universe_ptr->parameters.g_mean_maternal_immunity,
    Rcpp::Named("g_sum_maternal_immunity")=universe_ptr->parameters.g_sum_maternal_immunity,
    Rcpp::Named("g_total_mums")=universe_ptr->parameters.g_total_mums,
    Rcpp::Named("g_N")=universe_ptr->parameters.g_N,
    Rcpp::Named("g_max_age")=universe_ptr->parameters.g_max_age,
    Rcpp::Named("g_average_age")=universe_ptr->parameters.g_average_age,
    Rcpp::Named("g_EIR")=universe_ptr->parameters.g_EIR,
    Rcpp::Named("g_a0")=universe_ptr->parameters.g_a0,
    Rcpp::Named("g_rho")=universe_ptr->parameters.g_rho,
    Rcpp::Named("g_zeta_meanlog")=universe_ptr->parameters.g_zeta_meanlog,
    Rcpp::Named("g_zeta_sdlog")=universe_ptr->parameters.g_zeta_sdlog,
    Rcpp::Named("g_ft")=universe_ptr->parameters.g_ft,
    Rcpp::Named("g_mu0")=universe_ptr->parameters.g_mu0,
    Rcpp::Named("g_mean_mosquito_age")=universe_ptr->parameters.g_mean_mosquito_age,
    Rcpp::Named("g_beta_gradient")=universe_ptr->parameters.g_beta_gradient,
    Rcpp::Named("g_beta_intercept")=universe_ptr->parameters.g_beta_intercept,
    Rcpp::Named("g_ak")=universe_ptr->parameters.g_ak,
    Rcpp::Named("g_Q0")=universe_ptr->parameters.g_Q0,
    Rcpp::Named("g_mosquito_deficit")=universe_ptr->parameters.g_mosquito_deficit,
    Rcpp::Named("g_scourge_today")=universe_ptr->parameters.g_scourge_today,
    Rcpp::Named("g_mean_mv")=universe_ptr->parameters.g_mean_mv,
    Rcpp::Named("g_mosquito_next_biting_day_vector")=universe_ptr->parameters.g_mosquito_next_biting_day_vector,
    Rcpp::Named("g_mosquito_biting_counter")=universe_ptr->parameters.g_mosquito_biting_counter,
    Rcpp::Named("g_max_mosquito_biting_counter")=universe_ptr->parameters.g_max_mosquito_biting_counter,
    Rcpp::Named("g_delay_mos")=universe_ptr->parameters.g_delay_mos,
    Rcpp::Named("g_delay_gam")=universe_ptr->parameters.g_delay_gam,
    Rcpp::Named("g_dur_E")=universe_ptr->parameters.g_dur_E,
    Rcpp::Named("g_dur_T")=universe_ptr->parameters.g_dur_T,
    Rcpp::Named("g_dur_D")=universe_ptr->parameters.g_dur_D,
    Rcpp::Named("g_dur_U")=universe_ptr->parameters.g_dur_U,
    Rcpp::Named("g_dur_P")=universe_ptr->parameters.g_dur_P,
    Rcpp::Named("g_dur_A")=universe_ptr->parameters.g_dur_A,
    Rcpp::Named("g_dur_AU")=universe_ptr->parameters.g_dur_AU,
    Rcpp::Named("g_d1")=universe_ptr->parameters.g_d1,
    Rcpp::Named("g_dID")=universe_ptr->parameters.g_dID,
    Rcpp::Named("g_ID0")=universe_ptr->parameters.g_ID0,
    Rcpp::Named("g_kD")=universe_ptr->parameters.g_kD,
    Rcpp::Named("g_uD")=universe_ptr->parameters.g_uD,
    Rcpp::Named("g_aD")=universe_ptr->parameters.g_aD,
    Rcpp::Named("g_fD0")=universe_ptr->parameters.g_fD0,
    Rcpp::Named("g_gD")=universe_ptr->parameters.g_gD,
    Rcpp::Named("g_alphaU")=universe_ptr->parameters.g_alphaU,
    Rcpp::Named("g_b0")=universe_ptr->parameters.g_b0,
    Rcpp::Named("g_b1")=universe_ptr->parameters.g_b1,
    Rcpp::Named("g_dB")=universe_ptr->parameters.g_dB,
    Rcpp::Named("g_IB0")=universe_ptr->parameters.g_IB0,
    Rcpp::Named("g_kB")=universe_ptr->parameters.g_kB,
    Rcpp::Named("g_uB")=universe_ptr->parameters.g_uB,
    Rcpp::Named("g_phi0")=universe_ptr->parameters.g_phi0,
    Rcpp::Named("g_phi1")=universe_ptr->parameters.g_phi1,
    Rcpp::Named("g_dCA")=universe_ptr->parameters.g_dCA,
    Rcpp::Named("g_IC0")=universe_ptr->parameters.g_IC0,
    Rcpp::Named("g_kC")=universe_ptr->parameters.g_kC,
    Rcpp::Named("g_uCA")=universe_ptr->parameters.g_uCA,
    Rcpp::Named("g_PM")=universe_ptr->parameters.g_PM,
    Rcpp::Named("g_dCM")=universe_ptr->parameters.g_dCM,
    Rcpp::Named("g_gamma1")=universe_ptr->parameters.g_gamma1,
    Rcpp::Named("g_cD")=universe_ptr->parameters.g_cD,
    Rcpp::Named("g_cT")=universe_ptr->parameters.g_cT,
    Rcpp::Named("g_cU")=universe_ptr->parameters.g_cU,
    Rcpp::Named("g_total_human_infections")=universe_ptr->parameters.g_total_human_infections,
    Rcpp::Named("g_total_mosquito_infections")=universe_ptr->parameters.g_total_mosquito_infections,
    Rcpp::Named("g_percentage_imported_human_infections")=universe_ptr->parameters.g_percentage_imported_human_infections,
    Rcpp::Named("g_percentage_imported_mosquito_infections")=universe_ptr->parameters.g_percentage_imported_mosquito_infections,
    Rcpp::Named("g_cotransmission_frequencies")=universe_ptr->parameters.g_cotransmission_frequencies,
    Rcpp::Named("g_cotransmission_frequencies_counter")=universe_ptr->parameters.g_cotransmission_frequencies_counter,
    Rcpp::Named("g_cotransmission_frequencies_size")=universe_ptr->parameters.g_cotransmission_frequencies_size,
    Rcpp::Named("g_oocyst_frequencies")=universe_ptr->parameters.g_oocyst_frequencies,
    Rcpp::Named("g_oocyst_frequencies_counter")=universe_ptr->parameters.g_oocyst_frequencies_counter,
    Rcpp::Named("g_oocyst_frequencies_size")=universe_ptr->parameters.g_oocyst_frequencies_size,
    // genetics
    Rcpp::Named("g_identity_id")=Parameters::g_identity_id,
    Rcpp::Named("g_num_loci")=Parameters::g_num_loci,
    Rcpp::Named("g_ibd_length")=Parameters::g_ibd_length,
    Rcpp::Named("g_barcode_length")=Parameters::g_barcode_length,
    Rcpp::Named("g_plaf")=Parameters::g_plaf,
    Rcpp::Named("g_prob_crossover")=Parameters::g_prob_crossover,
    Rcpp::Named("g_barcode_type")=static_cast<unsigned int>(Parameters::g_barcode_type),
    Rcpp::Named("g_island_imports_plaf_linked_flag")=Parameters::g_island_imports_plaf_linked_flag,
    // spatial
    Rcpp::Named("g_spatial_type")=static_cast<unsigned int>(universe_ptr->parameters.g_spatial_type),
    Rcpp::Named("g_spatial_imported_cotransmission_frequencies")=universe_ptr->parameters.g_spatial_imported_cotransmission_frequencies,
    Rcpp::Named("g_spatial_imported_oocyst_frequencies")=universe_ptr->parameters.g_spatial_imported_oocyst_frequencies,
    Rcpp::Named("g_spatial_exported_cotransmission_frequencies")=universe_ptr->parameters.g_spatial_exported_cotransmission_frequencies,
    Rcpp::Named("g_spatial_exported_oocyst_frequencies")=universe_ptr->parameters.g_spatial_exported_oocyst_frequencies,
    Rcpp::Named("g_spatial_total_exported_barcodes")=universe_ptr->parameters.g_spatial_total_exported_barcodes,
    Rcpp::Named("g_spatial_total_exported_oocysts")=universe_ptr->parameters.g_spatial_total_exported_oocysts,
    Rcpp::Named("g_spatial_total_imported_human_infections")=universe_ptr->parameters.g_spatial_total_imported_human_infections,
    Rcpp::Named("g_spatial_total_imported_mosquito_infections")=universe_ptr->parameters.g_spatial_total_imported_mosquito_infections,
    Rcpp::Named("g_spatial_imported_human_infection_counter")=universe_ptr->parameters.g_spatial_imported_human_infection_counter,
    Rcpp::Named("g_spatial_imported_mosquito_infection_counter")=universe_ptr->parameters.g_spatial_imported_mosquito_infection_counter,
    Rcpp::Named("g_spatial_exported_barcode_counter")=universe_ptr->parameters.g_spatial_exported_barcode_counter,
    Rcpp::Named("g_spatial_exported_oocyst_counter")=universe_ptr->parameters.g_spatial_exported_oocyst_counter,
    Rcpp::Named("g_spatial_exported_barcodes")=spatial_exported_barcodes,
    Rcpp::Named("g_spatial_imported_barcodes")=spatial_imported_barcodes,
    Rcpp::Named("g_spatial_exported_oocysts")=spatial_exported_oocysts,
    Rcpp::Named("g_spatial_imported_oocysts")=spatial_imported_oocysts,
    // Drugs and resistance
    Rcpp::Named("g_resistance_flag")=universe_ptr->parameters.g_resistance_flag,
    Rcpp::Named("g_absolute_fitness_cost_flag")=universe_ptr->parameters.g_absolute_fitness_cost_flag,
    Rcpp::Named("g_number_of_resistance_loci")=universe_ptr->parameters.g_number_of_resistance_loci,
    Rcpp::Named("g_resistance_loci")=universe_ptr->parameters.g_resistance_loci,
    Rcpp::Named("g_artemisinin_loci")=universe_ptr->parameters.g_artemisinin_loci,
    Rcpp::Named("g_cost_of_resistance")=universe_ptr->parameters.g_cost_of_resistance,
    Rcpp::Named("g_drugs")=drugs,
    Rcpp::Named("g_mft_flag")=universe_ptr->parameters.g_mft_flag,
    Rcpp::Named("g_partner_drug_ratios")=universe_ptr->parameters.g_partner_drug_ratios,
    Rcpp::Named("g_drug_choice")=universe_ptr->parameters.g_drug_choice,
    Rcpp::Named("g_number_of_drugs")=universe_ptr->parameters.g_number_of_drugs,
    Rcpp::Named("g_dur_SPC")=universe_ptr->parameters.g_dur_SPC,
    // mutation
    Rcpp::Named("g_mutation_flag")=universe_ptr->parameters.g_mutation_flag,
    Rcpp::Named("g_mutation_rate")=universe_ptr->parameters.g_mutation_rate,
    Rcpp::Named("g_mutation_treated_modifier")=universe_ptr->parameters.g_mutation_treated_modifier,
    Rcpp::Named("g_mutations_today")=universe_ptr->parameters.g_mutations_today,
    Rcpp::Named("g_mutation_pos_allocator")=universe_ptr->parameters.g_mutation_pos_allocator,
    // vector adaptation
    Rcpp::Named("g_vector_adaptation_flag")=universe_ptr->parameters.g_vector_adaptation_flag,
    Rcpp::Named("g_vector_adaptation_loci")=universe_ptr->parameters.g_vector_adaptation_loci,
    Rcpp::Named("g_local_oocyst_advantage")=universe_ptr->parameters.g_local_oocyst_advantage,
    Rcpp::Named("g_gametocyte_sterilisation_flag")=universe_ptr->parameters.g_gametocyte_sterilisation_flag,
    Rcpp::Named("g_gametocyte_sterilisation")=universe_ptr->parameters.g_gametocyte_sterilisation,
    Rcpp::Named("g_oocyst_reduction_by_artemisinin")=universe_ptr->parameters.g_oocyst_reduction_by_artemisinin,
    // nmf
    Rcpp::Named("g_nmf_flag")=universe_ptr->parameters.g_nmf_flag,
    Rcpp::Named("g_mean_nmf_frequency")=universe_ptr->parameters.g_mean_nmf_frequency,
    Rcpp::Named("g_nmf_age_brackets")=universe_ptr->parameters.g_nmf_age_brackets,
    Rcpp::Named("g_prob_of_testing_nmf")=universe_ptr->parameters.g_prob_of_testing_nmf
  );
  
  Rcpp::List parameter_housekeeping_List = Rcpp::List::create(
    // housekeeping
    Rcpp::Named("g_h_quiet_print")=universe_ptr->parameters.g_h_quiet_print,
    Rcpp::Named("g_h_quiet_test_print")=universe_ptr->parameters.g_h_quiet_test_print
  );
  
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // END: MODEL STATE GET
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
  auto duration = chrono::duration_cast<std::chrono::seconds>(t1 - t0).count();
   rcpp_out(universe_ptr->parameters.g_h_quiet_print, "Time elapsed in fetching state: " + std::to_string(duration) + " seconds\n");
  
  // Return Named List with population and parameters
  return Rcpp::List::create(
    Rcpp::Named("population_List") = population_List, 
    Rcpp::Named("populations_event_and_strains_List") = populations_event_and_strains_List,
    Rcpp::Named("scourge_List") = scourge_List,
    Rcpp::Named("parameters_List")=parameters_List,
    Rcpp::Named("parameter_housekeeping_List")=parameter_housekeeping_List
    );
  
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // fini
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
}


