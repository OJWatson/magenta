//
//  MAGENTA
//  main_get.cpp
//
//  Created: OJ Watson on 06/12/2015
//
//  Distributed under the MIT software licence
//
//  Further details (if any) of this set of functions can be found in the corresponding header file.
//
// ---------------------------------------------------------------------------

#include <RcppArmadillo.h>
#include "stdafx.h"
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
// [[Rcpp::export]]
Rcpp::List Simulation_Get_cpp(Rcpp::List paramList)
{
  
  // prove that C++ code is being run
  Rcpp::Rcout << "Rcpp function is working!\n";
  
  // start timer
  chrono::high_resolution_clock::time_point t0 = std::chrono::high_resolution_clock::now();
  
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // START: R -> C++ CONVERSIONS
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  // Create universe pointer from paramList statePtr
  Rcpp::XPtr<Universe> universe_ptr = Rcpp::as<Rcpp::XPtr<Universe> > (paramList["statePtr"]);
  
  // Initialise variables from the statePtr provided that are needing to be saved
  std::vector<double> Zeta = universe_ptr->zeta_vector;
  
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // END: R -> C++ CONVERSIONS
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  Rcpp::Rcout << "Pointer unpacking working!\n";
  
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
  std::vector<double> ID(universe_ptr->parameters.g_N);
  std::vector<int> IB_last_boost_time(universe_ptr->parameters.g_N);
  std::vector<int> ICA_last_boost_time(universe_ptr->parameters.g_N);
  std::vector<int> ID_last_boost_time(universe_ptr->parameters.g_N);
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
  std::vector<std::vector<int> > Infection_time_realisation_vectors(universe_ptr->parameters.g_N);
  std::vector<std::vector<int> > Infection_state_realisation_vectors(universe_ptr->parameters.g_N);
  std::vector<std::vector<std::vector<bool> > > Infection_barcode_realisation_vectors(universe_ptr->parameters.g_N);
  std::vector<std::vector<std::vector<bool> > > Strain_barcode_vectors(universe_ptr->parameters.g_N);
  std::vector<std::vector<int> > Strain_infection_state_vectors(universe_ptr->parameters.g_N);
  std::vector<std::vector<int> > Strain_day_of_infection_state_change_vectors(universe_ptr->parameters.g_N);
  
  Rcpp::Rcout << "Vector initialisation working!\n";
  
  // Temporary necessities for casting vectors for pending states
  std::vector<Person::InfectionStatus> temp_infection_state_realisation_vector{};
  unsigned int temp_status_iterator = 0;
  
  // Temporary necessities for pending barcodes
  std::vector<barcode_t> temp_infection_barcode_realisation_vector = {};
  std::vector<bool> temp_barcode_bool_vector{};
  temp_barcode_bool_vector.reserve(barcode_length);
  int temp_barcode_iterator = 0;
  
  // Temporary necessities for strains
  Strain temp_strain;
  barcode_t temp_barcode;
  int temp_strain_iterator = 0;

  Rcpp::Rcout << "Preloop working!\n";

  
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
    
    // Strain Numbers
    Number_of_Strains[element] = universe_ptr->population[element].get_m_number_of_strains();
    
    // Temp strain to next change state
    Day_of_next_strain_state_change[element] = universe_ptr->population[element].get_m_day_of_next_strain_state_change();
    Day_of_next_event[element] = universe_ptr->population[element].get_m_day_of_next_event();
    
    // Realised Infections
    Number_of_Realised_Infections[element] = universe_ptr->population[element].get_m_number_of_realised_infections();
    
    // Pending Infection time vector
    Infection_time_realisation_vectors[element] = universe_ptr->population[element].get_m_infection_time_realisation_vector();
    
    //Rcpp::Rcout << "Prepending working!\n";
    
    // Pending Infection barcode and state vector
    // ---------------------------------------
    temp_infection_barcode_realisation_vector = universe_ptr->population[element].get_m_infection_barcode_realisation_vector();
    Infection_barcode_realisation_vectors[element].reserve(temp_infection_barcode_realisation_vector.size());
   
    temp_infection_state_realisation_vector = universe_ptr->population[element].get_m_infection_state_realisation_vector();
    Infection_state_realisation_vectors[element].reserve(temp_infection_state_realisation_vector.size());
    
    for(temp_status_iterator = 0 ; temp_status_iterator < temp_infection_state_realisation_vector.size() ; temp_status_iterator++)
    {
      Infection_state_realisation_vectors[element].push_back(temp_infection_state_realisation_vector[temp_status_iterator]);
      // fetch barcode and turn into vector<bool>
      for(temp_barcode_iterator = 0; temp_barcode_iterator < barcode_length ; temp_barcode_iterator++ )
      {
        temp_barcode_bool_vector.push_back(temp_infection_barcode_realisation_vector[temp_status_iterator][temp_barcode_iterator]);
      }
      Infection_barcode_realisation_vectors[element].push_back(temp_barcode_bool_vector);
      temp_barcode_bool_vector.clear();

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
    
    // Loop through each strain converting into barcodes, states and days of infection state changes
    for(temp_strain_iterator = 0 ; temp_strain_iterator < Number_of_Strains[element] ; temp_strain_iterator++)
    {
      
      temp_strain = universe_ptr->population[element].get_m_person_strain_x(temp_strain_iterator);
      
      Strain_infection_state_vectors[element].push_back(temp_strain.get_m_strain_infection_status());
      Strain_day_of_infection_state_change_vectors[element].push_back(temp_strain.get_m_day_of_strain_infection_status_change());
      // fetch barcode and turn into vector<bool>
      temp_barcode = temp_strain.get_m_barcode();
      for(temp_barcode_iterator = 0; temp_barcode_iterator < barcode_length ; temp_barcode_iterator++ )
      {
        temp_barcode_bool_vector.push_back(temp_barcode[temp_barcode_iterator]);
      }
      Strain_barcode_vectors[element].push_back(temp_barcode_bool_vector);
      temp_barcode_bool_vector.clear();
  
    }
    
  }
    
    // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    // END: HUMAN POPULATION
    // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    
    // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    // START: MOSQUITO POPULATION
    // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    
    // Create model state storage
    unsigned short int scourge_size = universe_ptr->scourge.size();
  
    std::vector<int> Mosquito_Infection_States(scourge_size);
    std::vector<unsigned short int> Mosquito_Day_of_next_blood_meal(scourge_size);		
    std::vector<unsigned short int> Mosquito_Day_of_death(scourge_size);
    std::vector<unsigned short int> Mosquito_Number_of_ruptured_oocysts (scourge_size);
    std::vector<std::vector<unsigned short int> > Mosquito_Oocyst_rupture_time_vectors(scourge_size);
    std::vector<std::vector<std::vector<bool> > > Mosquito_Oocyst_barcode_male_vectors(scourge_size);
    std::vector<std::vector<std::vector<bool> > > Mosquito_Oocyst_barcode_female_vectors(scourge_size);
    
    // Temporary necessities for pending barcodes
    std::vector<barcode_t> temp_male_barcode_realisation_vector = {};
    std::vector<barcode_t> temp_female_barcode_realisation_vector = {};

    std::vector<bool> temp_barcode_male_bool_vector{};
    std::vector<bool> temp_barcode_female_bool_vector{};
    temp_barcode_female_bool_vector.reserve(barcode_length);
    temp_barcode_male_bool_vector.reserve(barcode_length);
    
    Rcpp::Rcout << "Premosquito loop working!\n";

    for (unsigned int element = 0; element < scourge_size ; element++) 
    {
      
      // Infection States
      Mosquito_Infection_States[element] = static_cast<int>(universe_ptr->scourge[element].get_m_mosquito_infection_state());

      // Day Changes
      Mosquito_Day_of_death[element] = universe_ptr->scourge[element].get_m_day_of_death();

      Mosquito_Day_of_next_blood_meal[element] = universe_ptr->scourge[element].get_m_day_of_next_blood_meal();

      // Ruptured oocyst Numbers 
      Mosquito_Number_of_ruptured_oocysts[element] = static_cast<unsigned short int>(universe_ptr->scourge[element].get_m_ruptured_oocyst_count());

      // Pending Infection time vector
      Mosquito_Oocyst_rupture_time_vectors[element] = universe_ptr->scourge[element].get_m_oocyst_rupture_time_vector();

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
        for(temp_barcode_iterator = 0; temp_barcode_iterator < barcode_length ; temp_barcode_iterator++ )
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
    
  
  
  Rcpp::Rcout << "Postloop working!\n";
  
  // Create Rcpp population list
  Rcpp::List population_List = Rcpp::List::create(
    Rcpp::Named("Infection_States")=Infection_States,
    Rcpp::Named("Zetas")=Zeta,
    Rcpp::Named("Ages")=Ages,
    Rcpp::Named("IB")=IB,
    Rcpp::Named("ICA")=ICA,
    Rcpp::Named("ICM")=ICM,
    Rcpp::Named("ID")=ID,
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
    Rcpp::Named("Number_of_Realised_Infections")=Number_of_Realised_Infections,
    Rcpp::Named("Infection_time_realisation_vectors")=Infection_time_realisation_vectors,
    Rcpp::Named("Infection_state_realisation_vectors")=Infection_state_realisation_vectors,
    Rcpp::Named("Infection_barcode_realisation_vectors")=Infection_barcode_realisation_vectors,
    Rcpp::Named("Strain_infection_state_vectors")=Strain_infection_state_vectors,
    Rcpp::Named("Strain_day_of_infection_state_change_vectors")=Strain_day_of_infection_state_change_vectors,
    Rcpp::Named("Strain_barcode_vectors")=Strain_barcode_vectors
  );
  
  // Create Mosquito population list
  Rcpp::List scourge_List = Rcpp::List::create(
    Rcpp::Named("Mosquito_Infection_States")=Mosquito_Infection_States,
    Rcpp::Named("Mosquito_Day_of_next_blood_meal")=Mosquito_Day_of_next_blood_meal,
    Rcpp::Named("Mosquito_Day_of_death")=Mosquito_Day_of_death,
    Rcpp::Named("Mosquito_Number_of_ruptured_oocysts")=Mosquito_Number_of_ruptured_oocysts,
    Rcpp::Named("Mosquito_Oocyst_rupture_time_vectors")=Mosquito_Oocyst_rupture_time_vectors,
    Rcpp::Named("Mosquito_Oocyst_barcode_male_vectors")=Mosquito_Oocyst_barcode_male_vectors,
    Rcpp::Named("Mosquito_Oocyst_barcode_female_vectors")=Mosquito_Oocyst_barcode_female_vectors,
    Rcpp::Named("Scourge_size")=universe_ptr->scourge.size()
  );
  
  // Create Rcpp Parameters list
  Rcpp::List parameters_List = Rcpp::List::create(
    Rcpp::Named("g_current_time")=universe_ptr->parameters.g_current_time,
    Rcpp::Named("g_mean_maternal_immunity")=universe_ptr->parameters.g_mean_maternal_immunity,
    Rcpp::Named("g_sum_maternal_immunity")=universe_ptr->parameters.g_sum_maternal_immunity,
    Rcpp::Named("g_total_mums")=universe_ptr->parameters.g_total_mums,
    Rcpp::Named("g_N")=universe_ptr->parameters.g_N
  );
  
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // END: MODEL STATE GET
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
  auto duration = chrono::duration_cast<std::chrono::seconds>(t1 - t0).count();
  std::cout << "Time elapsed in fetching state: " << duration << " seconds" << std::endl;
  
  // Return Named List with population and parameters
  return Rcpp::List::create(
    Rcpp::Named("population_List") = population_List, 
    Rcpp::Named("populations_event_and_strains_List") = populations_event_and_strains_List,
    Rcpp::Named("scourge_List") = scourge_List,
    Rcpp::Named("parameters_List")=parameters_List);
  
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // fini
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
}


