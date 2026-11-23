//
//  magenta
//  main_update.cpp
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
#include "util.h"
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
//' Continues simulation forward for as long as specified in param_list
//'
//' @param param_list parameter list generated with \code{Param_List_Simulation_Update_Create}
//' @return list with ptr to model state and loggers describing the current model state
//' @export
// [[Rcpp::export]]
Rcpp::List Simulation_Update_cpp(Rcpp::List param_list)
{
  
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // START: R -> C++ CONVERSIONS
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  // Create universe pointer from param_list statePtr
  Rcpp::XPtr<Universe> u_ptr = Rcpp::as<Rcpp::XPtr<Universe> > (param_list["statePtr"]);
  
  // prove that C++ code is being run
  rcpp_out(u_ptr->parameters.g_h_quiet_print, "Rcpp function is working!\n");
  
  // Initialise all the universal variables from the statePtr provided
  u_ptr->parameters.g_years = Rcpp::as<double>(param_list["years"]);
  u_ptr->parameters.g_ft = Rcpp::as<double>(param_list["ft"]);
  Rcpp::List spatial_list = param_list["spatial_list"];
  Rcpp::List drug_list = param_list["drug_list"];
  Rcpp::List barcode_list  = param_list["barcode_list"];
  
  // Spatial updates 
  // Metapopulation not fully implemented yet
  if(u_ptr->parameters.g_spatial_type == Parameters::METAPOPULATION)
  {
    
    // convert R vector
    std::vector<std::vector<bool> > x = (Rcpp::as<std::vector<std::vector<bool> > >(param_list["imported_barcodes"]));
    
    // prep import barcode parameters and fill accordingly 
    u_ptr->parameters.g_spatial_imported_barcodes.reserve(x.size());
    
    // generate temp barcode iterator
    unsigned int temp_barcode_iterator = 0;
    
    // loop through imported barcodes and set
    for (unsigned int i = 0; i < x.size(); i++) 
    {
      // fetch vector<bool> and turn into barcode
      for (temp_barcode_iterator = 0; temp_barcode_iterator < Parameters::g_barcode_length; temp_barcode_iterator++)
      {
        Strain::temp_barcode[temp_barcode_iterator] = x[i][temp_barcode_iterator];
      }
      u_ptr->parameters.g_spatial_imported_barcodes[i] = Strain::temp_barcode;
    }
    
    // 
    
  } 
  else if (u_ptr->parameters.g_spatial_type == Parameters::ISLAND) 
  {
    u_ptr->parameters.g_percentage_imported_human_infections = Rcpp::as<double>(spatial_list["imported_cotransmissions_events"]);
    u_ptr->parameters.g_percentage_imported_mosquito_infections = Rcpp::as<double>(spatial_list["imported_oocyst_events"]);
    u_ptr->parameters.g_spatial_total_imported_human_infections = u_ptr->parameters.g_percentage_imported_human_infections * static_cast<double>(u_ptr->parameters.g_total_human_infections);
    u_ptr->parameters.g_spatial_total_imported_mosquito_infections = u_ptr->parameters.g_percentage_imported_mosquito_infections * static_cast<double>(u_ptr->parameters.g_total_mosquito_infections);
    u_ptr->parameters.g_spatial_total_imported_human_infections = rpoisson1(u_ptr->parameters.g_spatial_total_imported_human_infections);
    u_ptr->parameters.g_spatial_total_imported_mosquito_infections = rpoisson1(u_ptr->parameters.g_spatial_total_imported_mosquito_infections);
  }
  
  // Other updates
  u_ptr->parameters.g_plaf = Rcpp::as<std::vector<double> >(spatial_list["plaf"]);
  u_ptr->parameters.g_mutation_flag = Rcpp::as<bool>(barcode_list["mutation_flag"]);
  
  // Resistance updates
  if (drug_list["resistance_flag"]) {
    u_ptr->parameters.g_resistance_flag = Rcpp::as<bool>(drug_list["resistance_flag"]);
    u_ptr->parameters.g_drug_choice = Rcpp::as<int>(drug_list["drug_choice"]);
    u_ptr->parameters.g_partner_drug_ratios = Rcpp::as<std::vector<double> >(drug_list["partner_drug_ratios"]);
  }
  
  u_ptr->parameters.g_cotransmission_frequencies = Rcpp::as<std::vector<int> >(spatial_list["cotransmission_freq_vector"]);
  u_ptr->parameters.g_cotransmission_frequencies_size = u_ptr->parameters.g_cotransmission_frequencies.size();
  u_ptr->parameters.g_oocyst_frequencies = Rcpp::as<std::vector<int> >(spatial_list["oocyst_freq_vector"]);
  u_ptr->parameters.g_oocyst_frequencies_size = u_ptr->parameters.g_oocyst_frequencies.size();
  
  // Initialise the mosquito intervention parameter vectors
  // These vectors detail how interventions have had an effect on mosquito behaviour for the update length considered
  std::vector<double> mosquito_death_rates = Rcpp::as<vector<double> >(param_list["mu_vec"]);
  std::vector<double> mosquito_biting_rates = Rcpp::as<vector<double> >(param_list["fv_vec"]);
  
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // END: R -> C++ CONVERSIONS
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  rcpp_out(u_ptr->parameters.g_h_quiet_print, "Pointer unpacking working!\n");
  double Scount = 0;
  for(auto &element : u_ptr->population){
    if(element.get_m_infection_state()==Person::SUSCEPTIBLE){
      Scount++;
    }
  }
  rcpp_out(u_ptr->parameters.g_h_quiet_print, "Starting susceptible population: " + std::to_string(Scount/u_ptr->parameters.g_N) + "\n");
  
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // START: PRE LOOP
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

  // End of simulation time
  int g_end_time = u_ptr->parameters.g_current_time + (u_ptr->parameters.g_years * 365);
  
  // Preallocations;
  unsigned int num_bites = 0;
  double psi_sum = 0;
  unsigned int scourge_size = u_ptr->scourge.size();
  unsigned int temp_deficit = 0;
  int intervention_counter = 0;
  double biting_log_q = 0;
  
  // status eq for logging and other logging variables
  unsigned int not_treated = 0;
  std::vector<unsigned int> succesful_treatments_by_drug(u_ptr->parameters.g_number_of_drugs,0);
  std::vector<unsigned int> unsuccesful_treatments_by_drug(u_ptr->parameters.g_number_of_drugs,0);
  std::vector<double> status_eq = { 0,0,0,0,0,0 };
  std::vector<double> mosq_status_eq = { 0,0,0 };
  unsigned int log_counter = 0;
  double total_incidence = 0;
  double total_incidence_05 = 0;
  int daily_incidence_return = 0;
  double total_severe = 0;
  double total_severe_05 = 0;
  int daily_bite_counters = 0;
  int daily_infectious_bite_counters = 0;
  std::vector<unsigned int> daily_mutations_per_loci(u_ptr->parameters.g_num_loci,0);
  
  // Maternal 
  double mean_psi = 0;
  double pi_sum = 0;
  
  // Bites
  std::vector<int> mosquito_biting_queue;
  mosquito_biting_queue.reserve(scourge_size);
  std::vector<int> bite_storage_queue;
  bite_storage_queue.reserve(scourge_size);
  IntegerVector bite_allocation(u_ptr->parameters.g_N);
  // biting allocation
  std::vector<double> bite_randoms(scourge_size / 2);
  std::generate(bite_randoms.begin(), bite_randoms.end(), runif0_1);
  
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // START SIMULATION LOOP
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  for ( ; u_ptr->parameters.g_current_time < g_end_time ; u_ptr->parameters.g_current_time++, intervention_counter++)
  {
    
    // update the calendar day
    u_ptr->parameters.g_calendar_day++;
    if (u_ptr->parameters.g_calendar_day == 366 ) u_ptr->parameters.g_calendar_day = 1;
    
    
    // Counter print
    if (u_ptr->parameters.g_current_time % 100 == 0) 
    { 
      rcpp_out(u_ptr->parameters.g_h_quiet_print, std::to_string(u_ptr->parameters.g_current_time) + " days" + "\n"); 
    }
    
    // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    // DAILY UPDATING AND EVENT HANDLING
    // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    
    // First calculate the mean maternal immunity from yesterday
    // --------------------------------------------------------------------------------------------------------------------------------------------------
    
    // Calculate yesterday's mean maternal immunity
    u_ptr->parameters.g_mean_maternal_immunity = u_ptr->parameters.g_sum_maternal_immunity / u_ptr->parameters.g_total_mums;
    u_ptr->parameters.g_mean_maternal_severe_immunity = u_ptr->parameters.g_sum_maternal_severe_immunity / u_ptr->parameters.g_total_mums;
    
    // Reset maternal immunity sums
    u_ptr->parameters.g_sum_maternal_immunity = 0;
    u_ptr->parameters.g_sum_maternal_severe_immunity = 0;
    u_ptr->parameters.g_total_mums = 0;
    
    // Second reset counters
    // --------------------------------------------------------------------------------------------------------------------------------------------------
    
    // Reset import counters
    if (u_ptr->parameters.g_spatial_type == Parameters::ISLAND) 
    {
      u_ptr->parameters.g_spatial_total_imported_human_infections = u_ptr->parameters.g_percentage_imported_human_infections * u_ptr->parameters.g_total_human_infections;
      u_ptr->parameters.g_spatial_total_imported_mosquito_infections = u_ptr->parameters.g_percentage_imported_mosquito_infections * u_ptr->parameters.g_total_mosquito_infections;
    }
    u_ptr->parameters.g_spatial_imported_human_infection_counter = 0;
    u_ptr->parameters.g_spatial_imported_mosquito_infection_counter = 0;
    
    // Reset mutation counters
    if(u_ptr->parameters.g_mutation_flag){
        
      u_ptr->parameters.g_mutation_pos_allocator = 0;
      
      // mutation updates
      for(unsigned int l = 0; l < u_ptr->parameters.g_num_loci; l++){
        
        // if the mutation treated modifier is 1 then the same mutation rate is assumed regardless of the 
        // infection outcome, i.e. treated vs not treated. If this is true we can use pre-calculated mutations
        if(u_ptr->parameters.g_mutation_treated_modifier == 1.0) {
          
        u_ptr->parameters.g_mutations_today[l] = rpoisson1(u_ptr->parameters.g_mutation_rate[l] * u_ptr->parameters.g_total_human_infections);
        daily_mutations_per_loci[l] = daily_mutations_per_loci[l] + u_ptr->parameters.g_mutations_today[l];
        
        } else {
          u_ptr->parameters.g_mutations_today[l] = 0;
        }
      }
      
    }
    
    // Reset age dependent biting rate sum
    psi_sum = 0;
    
    // Reset total infections
    u_ptr->parameters.g_total_human_infections = 0;
    u_ptr->parameters.g_total_mosquito_infections = 0;
    
    // Third calculate the average biting rate and mosquito mortality for today given interventions
    // --------------------------------------------------------------------------------------------------------------------------------------------------
    
    // mortality changes due to interventions so we set it each day
    u_ptr->parameters.g_mean_mosquito_age = 1.0/mosquito_death_rates[intervention_counter];
    
    // frequency of biting also changes due to repel effects of interventions so generate next biting times.
    // These are geometric draws (days until next feed, at least 1) with mean 1/fv, i.e. each day a mosquito
    // feeds with probability fv, which is the daily step equivalent of the deterministic model's biting rate.
    // A fixed interval (e.g. always 3 days) meant the first bite after the 10 day EIP was always on day 12,
    // which gave fewer infectious bites per mosquito infection than the deterministic model
    biting_log_q = std::log(1.0 - std::min(mosquito_biting_rates[intervention_counter], 0.999999));
    for (auto &next_bite : u_ptr->parameters.g_mosquito_next_biting_day_vector) {
      next_bite = 1 + static_cast<int>(std::log(unif_rand()) / biting_log_q);
    }
    
    // EIP for mosquitoes infected today, drawn per bite as floor or ceiling of delay_mos + 0.5. On daily steps
    // a bite on the day of oocyst rupture is infectious, so an EIP of delay_mos would count a bite the
    // deterministic model only counts half of. Using delay_mos and delay_mos + 1 equally often gives the
    // deterministic model's expected number of infectious bites per mosquito infection
    for (auto &eip : u_ptr->parameters.g_mosquito_eip_vector) {
      eip = static_cast<int>(std::floor(u_ptr->parameters.g_delay_mos + 0.5));
      if (unif_rand() < (u_ptr->parameters.g_delay_mos + 0.5 - eip)) eip++;
    }
    
    
    // Loop through each person and mosquito and update
    // --------------------------------------------------------------------------------------------------------------------------------------------------
    
    // Human update loop
    // rcpp_out(u_ptr->parameters.g_h_quiet_print, "Human loop" + "\n"); 
    for (unsigned int human_update_i = 0; human_update_i < u_ptr->parameters.g_N; human_update_i++)
    {
      psi_sum += u_ptr->psi_vector[human_update_i] = u_ptr->population[human_update_i].update(u_ptr->parameters);
    }
    
    // Reset number of bites for each day and update each mosquito. Update returns whether the mosquito is biting today
    num_bites = 0;
    
    // Mosquito update loop
    
    // How many mosquitos ashould be on season for the day
    // Set deficit equal to how many were on season yesterday
    u_ptr->parameters.g_mosquito_deficit = u_ptr->parameters.g_scourge_today;
    
    // Calculate today's seasonal amount
    u_ptr->parameters.g_scourge_today = u_ptr->parameters.g_mean_mv * u_ptr->parameters.g_theta[ u_ptr->parameters.g_calendar_day -1];
    
    // Work out the difference
    u_ptr->parameters.g_mosquito_deficit -= u_ptr->parameters.g_scourge_today; // i.e. a positive deficit is a surplus of mosquitos
    
    temp_deficit = (u_ptr->parameters.g_mosquito_deficit > 0) ? (u_ptr->parameters.g_scourge_today + u_ptr->parameters.g_mosquito_deficit) : u_ptr->parameters.g_scourge_today;
    
    // loop through the mosquitos up to the number that should be active today
    // rcpp_out(u_ptr->parameters.g_h_quiet_print, "Mosquito loop" + "\n"); 
    for (unsigned int mosquito_update_i = 0; mosquito_update_i < temp_deficit; mosquito_update_i++)
    {
      // if the mosquito considered is more than needed then set to off season
      if (mosquito_update_i >= u_ptr->parameters.g_scourge_today) 
      {
        u_ptr->scourge[mosquito_update_i].set_m_mosquito_off_season(true);
        u_ptr->parameters.g_mosquito_deficit--;
      }
      else // otherwise update the mosquito
      {
        
        // if the mosquito is off season and we are considering it, i.e. we need it then set it equal to a random mosquito
        // and check if that mosquito is biting today by just looking at its next blood meal time as it will have already been 
        // updated.
        if(u_ptr->scourge[mosquito_update_i].get_m_mosquito_off_season()){
          u_ptr->scourge[mosquito_update_i] = u_ptr->scourge[runiform_int_1(0,mosquito_update_i-1)];
          u_ptr->parameters.g_mosquito_deficit++;
          if (u_ptr->scourge[mosquito_update_i].m_mosquito_biting_today) {
            mosquito_biting_queue.emplace_back(mosquito_update_i);
            num_bites++;
          }
        }
        else // otherwise just update it
        {
          if (u_ptr->scourge[mosquito_update_i].update(u_ptr->parameters))
          {
            mosquito_biting_queue.emplace_back(mosquito_update_i);
            num_bites++;
          }
        }
      }
    }
    
    // rcpp_out(u_ptr->parameters.g_h_quiet_print, "Shuffle bite queue" + "\n"); 
    // shuffle the bite queue otherwise you will introduce stepping-stone-esque genetic structuring
    shuffle_integer_vector(mosquito_biting_queue);
    
    // Adjust the number of bites to account for anthrophagy
    num_bites *= u_ptr->parameters.g_Q0;
    
    // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    // BITE HANDLING AND ALLOCATIONS
    // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    
    // First calculate mean age dependent biting heterogeneity (psi)
    // --------------------------------------------------------------------------------------------------------------------------------------------------
    
    // Calculate mean age dependent biting rate
    mean_psi = psi_sum / u_ptr->parameters.g_N;
    
    // Create normalised psi by dividing by the mean age dependent biting rate
    std::transform(u_ptr->psi_vector.begin(), u_ptr->psi_vector.end(), u_ptr->psi_vector.begin(),
                   std::bind(std::multiplies<double>(), 1 / mean_psi, std::placeholders::_1));
    
    // Create overall relative biting rate, pi, i.e. the product of individual biting heterogeneity and age dependent heterogeneity
    std::transform(u_ptr->psi_vector.begin(), u_ptr->psi_vector.end(),
                   u_ptr->zeta_vector.begin(), u_ptr->pi_vector.begin(),
                   std::multiplies<double>());
    
    // Caluclate total probability of being bitten within population
    pi_sum = std::accumulate(u_ptr->pi_vector.begin(), u_ptr->pi_vector.end(), 0.0);
    
    // Multinomial bite decision
    rmultinomN(num_bites, u_ptr->pi_vector, pi_sum, u_ptr->parameters.g_N, bite_storage_queue);
    
    // shuffle the biting queue of humans instead of mosquito queue as smaller and then saves spatial clustering from importation
    shuffle_integer_vector(bite_storage_queue);
    
    // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    // START: BITE ALLOCATIONS
    // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    
    // rcpp_out(u_ptr->parameters.g_h_quiet_print, "Bite allocation" + "\n"); 
    for (unsigned int num_bites_i = 0; num_bites_i < num_bites; num_bites_i++)
    {
      
      daily_bite_counters++;
      // allocate bite to human if mosquito is infected
      if (u_ptr->scourge[mosquito_biting_queue[num_bites_i]].get_m_mosquito_infection_state() == Mosquito::INFECTED) 
      {
        u_ptr->population[bite_storage_queue[num_bites_i]].allocate_bite(u_ptr->parameters, u_ptr->scourge[mosquito_biting_queue[num_bites_i]]);
        //if ( u_ptr->parameters.g_current_time > g_end_time - 7) 
        daily_infectious_bite_counters++;
      }
      
      // if human would cause infection to mosquito then handle the bite gametocytes
      if (u_ptr->population[bite_storage_queue[num_bites_i]].reciprocal_infection_boolean(u_ptr->parameters)) 
      {
        u_ptr->scourge[mosquito_biting_queue[num_bites_i]].handle_bite(u_ptr->parameters, u_ptr->population[bite_storage_queue[num_bites_i]]);
      }
      
    }
    
    // clear biting storage vectors
    bite_storage_queue.clear();
    mosquito_biting_queue.clear();
    
    
    // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    // START: SUMMARY LOGGING
    // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    
    // Mutation logging 
    // if the mutation treated modifier is not 1 then mutations need to be logged here as they are drawn at the infection rather than pre-computed
    if(u_ptr->parameters.g_mutation_treated_modifier != 1.0) {
    for(unsigned int l = 0; l < u_ptr->parameters.g_num_loci; l++){
        daily_mutations_per_loci[l] = daily_mutations_per_loci[l] + u_ptr->parameters.g_mutations_today[l];
      }
    }
    
    // Log the last period //
    
    if (u_ptr->parameters.g_current_time > g_end_time - (u_ptr->parameters.g_years * 365) - 1)
    {
      
      log_counter++;
      
      for (auto &element : u_ptr->population)
      {
        // Match infection state and schedule associated next state change
        switch (element.get_m_infection_state())
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
        
        // log incidence
        daily_incidence_return = element.log_daily_incidence(u_ptr->parameters);
        if(daily_incidence_return == 2)
        {
          total_incidence++;
          total_incidence_05++;
        } 
        if (daily_incidence_return == 1) 
        {
          total_incidence++;
        }
        
        // log expected severe cases from today's infections
        total_severe += element.get_m_expected_severe();
        if (element.get_m_person_age() <= (5 * 365))
        {
          total_severe_05 += element.get_m_expected_severe();
        }
        
        // log treatment outcomes
        if(element.get_m_treatment_outcome() == Person::SUCCESFULLY_TREATED){
          succesful_treatments_by_drug[element.get_m_drug_choice()]++;
        } 
        if (element.get_m_treatment_outcome() == Person::LPF){
          unsuccesful_treatments_by_drug[element.get_m_drug_choice()]++;
        }
        if (element.get_m_treatment_outcome() == Person::NOT_TREATED){
          not_treated++;
        }
      }
      
      for (auto &element : u_ptr->scourge)
      {
        // Match infection state and schedule associated next state change
        switch (element.get_m_mosquito_infection_state())
        {
        case Mosquito::SUSCEPTIBLE:
          mosq_status_eq[0]++;
          break;
        case Mosquito::EXPOSED :
          mosq_status_eq[1]++;
          break;
        case Mosquito::INFECTED:
          mosq_status_eq[2]++;
          break;
        default:
          assert(NULL && "Schedule Infection Status Change Error - mosquito's infection status not S, E, I");
        break;
        }
      }
    }
    
    // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    // END: SUMMARY LOGGING
    // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    
  };
  
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // SUMMARY LOGGING AVERAGING AND VARIABLE RETURN
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  
  // convert exproted barcodes to vector of vector of bool
  std::vector<SEXP>  Exported_Barcodes(u_ptr->parameters.g_spatial_total_exported_barcodes);
  if(u_ptr->parameters.g_spatial_type == Parameters::METAPOPULATION)
  {
    
    for(unsigned int temp_status_iterator = 0; temp_status_iterator < u_ptr->parameters.g_spatial_total_exported_barcodes ; temp_status_iterator++)
    {
      Exported_Barcodes[temp_status_iterator] = bitset_to_sexp(u_ptr->parameters.g_spatial_exported_barcodes[temp_status_iterator], u_ptr->parameters.g_barcode_length);
    }
    
  }
  
  // divide by population size and log counter and print to give overview
  rcpp_out(u_ptr->parameters.g_h_quiet_print, "S | D | A | U | T | P:\n");
  
  for (int element = 0; element < 6; element++) 
  {
    status_eq[element] /= (u_ptr->parameters.g_N * log_counter);
    rcpp_out(u_ptr->parameters.g_h_quiet_print, std::to_string(status_eq[element]) + " | ");
  }
  
  
  // Final infection states
  std::vector<int> Infection_States(u_ptr->parameters.g_N);
  std::vector<double> Ages(u_ptr->parameters.g_N);
  std::vector<double> IB(u_ptr->parameters.g_N);
  std::vector<double> ICA(u_ptr->parameters.g_N);
  std::vector<double> ICM(u_ptr->parameters.g_N);
  std::vector<double> ID(u_ptr->parameters.g_N);
  std::vector<double> IVA(u_ptr->parameters.g_N);
  std::vector<double> IVM(u_ptr->parameters.g_N);
  
  // loop through population and grabages and immunities for error checking
  for (unsigned int element = 0; element < u_ptr->parameters.g_N ; element++) 
  {
    
    // Ages and immunity 
    // TODO: Figure out the best way of standardising this logging 
    // Something like passing in a function name within the param_list which is the 
    // name for a logger written else where which then returns the Loggers obeject below
    Infection_States[element] = static_cast<int>(u_ptr->population[element].get_m_infection_state());
    u_ptr->population[element].update_immunities_to_today(u_ptr->parameters);
    Ages[element] = u_ptr->population[element].get_m_person_age();
    IB[element] = u_ptr->population[element].get_m_IB();
    ICA[element] = u_ptr->population[element].get_m_ICA();
    ICM[element] = u_ptr->population[element].get_m_ICM();
    ID[element] = u_ptr->population[element].get_m_ID();
    IVA[element] = u_ptr->population[element].get_m_IVA();
    IVM[element] = u_ptr->population[element].get_m_IVM();
    
  }
  
  rcpp_out(u_ptr->parameters.g_h_quiet_print, "\n Ages and Immunity Done \n");
  
  // Create Rcpp loggers list. Built with push_back as Rcpp::List::create takes at most 20 elements
  Rcpp::List Loggers;
  Loggers.push_back(log_counter, "Log_Counter");
  Loggers.push_back(total_incidence/log_counter, "Incidence");
  Loggers.push_back(total_incidence_05/log_counter, "Incidence_05");
  Loggers.push_back(total_severe/log_counter, "Severe");
  Loggers.push_back(total_severe_05/log_counter, "Severe_05");
  Loggers.push_back(Rcpp::List::create(
    Rcpp::Named("Successful_Treatments")=succesful_treatments_by_drug,
    Rcpp::Named("Unsuccesful_Treatments_LPF")=unsuccesful_treatments_by_drug,
    Rcpp::Named("Not_Treated")=not_treated,
    Rcpp::Named("daily_bite_counters")=daily_bite_counters,
    Rcpp::Named("daily_infectious_bite_counters")=daily_infectious_bite_counters
  ), "Treatments");
  Loggers.push_back(daily_mutations_per_loci, "daily_mutations_per_loci");
  Loggers.push_back(Infection_States, "InfectionStates");
  Loggers.push_back(Ages, "Ages");
  Loggers.push_back(IB, "IB");
  Loggers.push_back(ICA, "ICA");
  Loggers.push_back(ICM, "ICM");
  Loggers.push_back(ID, "ID");
  Loggers.push_back(IVA, "IVA");
  Loggers.push_back(IVM, "IVM");
  Loggers.push_back(Rcpp::List::create(
    Rcpp::Named("Mos_S")=mosq_status_eq[0]/log_counter,
    Rcpp::Named("Mos_E")=mosq_status_eq[1]/log_counter,
    Rcpp::Named("Mos_I")=mosq_status_eq[2]/log_counter
  ), "Mosquitoes");
  Loggers.push_back(u_ptr->parameters.g_current_time, "Time");
  Loggers.push_back(status_eq[0], "S");
  Loggers.push_back(status_eq[1], "D");
  Loggers.push_back(status_eq[2], "A");
  Loggers.push_back(status_eq[3], "U");
  Loggers.push_back(status_eq[4], "T");
  Loggers.push_back(status_eq[5], "P");
  
  rcpp_out(u_ptr->parameters.g_h_quiet_print, "\n Loggers Done \n");
  
  // Return Named List with pointer and loggers
  // If spatial also required then export the barcodes
  if(u_ptr->parameters.g_spatial_type == Parameters::METAPOPULATION)
  {
    return Rcpp::List::create(Rcpp::Named("Ptr") = u_ptr, Rcpp::Named("Loggers")=Loggers, Rcpp::Named("Exported_Barcodes")=Exported_Barcodes);
  } 
  else if(u_ptr->parameters.g_spatial_type == Parameters::ISLAND)
  {
    Rcpp::List Infections = Rcpp::List::create(Rcpp::Named("Humans")=u_ptr->parameters.g_total_human_infections,
                                               Rcpp::Named("Mosquitoes")=u_ptr->parameters.g_total_mosquito_infections);
    return Rcpp::List::create(Rcpp::Named("Ptr") = u_ptr, Rcpp::Named("Loggers")=Loggers, Rcpp::Named("Infections")=Infections);
  } 
  else 
  {
    return Rcpp::List::create(Rcpp::Named("Ptr") = u_ptr, Rcpp::Named("Loggers")=Loggers);
  }
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // fini
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
}


