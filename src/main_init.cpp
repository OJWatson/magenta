//
//  MAGENTA
//  main_init.cpp
//
//  Created: OJ Watson on 06/12/2015
//	Most recent edits: OJ on 10/03/2017
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
#include "person.h"
#include <cassert> // for error checking

#include <chrono>
#include <functional>
#include <numeric>  
#include <algorithm>
#include "util.h"

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
//' Creates initial model simulation using paramter list provided
//'
//' @param paramList parameter list generated with \code{Param_List_Simulation_Init_Create}
//' @return list with ptr to model state and loggers describing the current model state
//' @export
// [[Rcpp::export]]
Rcpp::List Simulation_Init_cpp(Rcpp::List paramList)
{
  
  // prove that C++ code is being run
  Rcpp::Rcout << "Rcpp function is working!\n";
  
  // Initialise parameters object
  Parameters parameters;
  
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // START: UNPACKING PARAMETERS
  //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  // Grab these R parameters first so that they are used in the initialisation
  parameters.g_N = Rcpp::as<unsigned int>(paramList["N"]);
  
  // Unpack the R equilibirum state parameter list object and barcode object
  Rcpp::List eqSS = paramList["eqSS"];
  Rcpp::List barcode_parms  = paramList["barcode_parms"];
  Rcpp::List spatial_list  = paramList["spatial_list"];
  
  // Un pack barcode parms
  parameters.g_num_loci = Rcpp::as<unsigned int>(barcode_parms["num_loci"]);
  parameters.g_ibd_length = Rcpp::as<unsigned int>(barcode_parms["ibd_length"]);
  parameters.g_barcode_length = static_cast<int>(parameters.g_num_loci * parameters.g_ibd_length);
  parameters.g_plaf = Rcpp::as<std::vector<double> >(barcode_parms["plaf"]);
  parameters.g_prob_crossover = Rcpp::as<std::vector<double> >(barcode_parms["prob_crossover"]);
  parameters.g_barcode_type = static_cast<Parameters::g_barcode_type_enum>(Rcpp::as<unsigned int>(barcode_parms["barcode_type"]));
  
  Strain::temp_barcode = boost::dynamic_bitset<>(Parameters::g_barcode_length);
  Strain::temp_identity_barcode = boost::dynamic_bitset<>(Parameters::g_ibd_length);
  Strain::temp_crossovers = boost::dynamic_bitset<>(Parameters::g_num_loci);
  
  // Grab seasonality and spatial
  parameters.g_spatial_type = static_cast<Parameters::g_spatial_type_enum>(Rcpp::as<unsigned int>(spatial_list["spatial_type"]));
  parameters.g_theta = Rcpp::as<vector<double> >(eqSS["theta"]);
  
  Rcpp::Rcout << "Running model of " << parameters.g_spatial_type << " spatial type.\n";

  // Mosquito steady state values at a population level
  double Sv = Rcpp::as<double>(eqSS["Sv"]);
  double Ev = Rcpp::as<double>(eqSS["Ev"]);
  double Iv = Rcpp::as<double>(eqSS["Iv"]);
   
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // START: INITIALISATION
  //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  // Add the human population
  std::vector<Person> population;
  population.reserve(parameters.g_N);
  
  // Initialise vectors for keeping biting related variables
  std::vector<double> psi_vector(parameters.g_N);
  std::vector<double> zeta_vector(parameters.g_N);
  std::vector<double> pi_vector(parameters.g_N);
  
  // Add the mosquito population, i.e. scourge
  std::vector<Mosquito> scourge;
  parameters.g_mean_mv = parameters.g_N * (Sv + Ev + Iv);
  scourge.reserve(static_cast<int>(parameters.g_N * (Sv + Ev + Iv) * *std::max_element(parameters.g_theta.begin(), parameters.g_theta.end())));
  parameters.g_scourge_today = parameters.g_mean_mv * parameters.g_theta[parameters.g_calendar_day];
  
  
  Rcpp::Rcout << "Mean mv = " << parameters.g_mean_mv << "!\n";
  Rcpp::Rcout << "Sv = " << Sv << "!\n";
  Rcpp::Rcout << "Ev = " << Ev << "!\n";
  Rcpp::Rcout << "Iv = " << Iv << "!\n";
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // END: INITIALISATION
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // START: R -> C++ CONVERSIONS
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  Rcpp::NumericMatrix Smat(Rcpp::as<NumericMatrix>(eqSS["Smat"]));
  Rcpp::NumericMatrix Dmat(Rcpp::as<NumericMatrix>(eqSS["Dmat"]));
  Rcpp::NumericMatrix Amat(Rcpp::as<NumericMatrix>(eqSS["Amat"]));
  Rcpp::NumericMatrix Umat(Rcpp::as<NumericMatrix>(eqSS["Umat"]));
  Rcpp::NumericMatrix Tmat(Rcpp::as<NumericMatrix>(eqSS["Tmat"]));
  Rcpp::NumericMatrix Pmat(Rcpp::as<NumericMatrix>(eqSS["Pmat"]));
  Rcpp::NumericMatrix IBmat(Rcpp::as<NumericMatrix>(eqSS["IBmat"]));
  Rcpp::NumericMatrix ICAmat(Rcpp::as<NumericMatrix>(eqSS["ICAmat"]));
  Rcpp::NumericMatrix ICMmat(Rcpp::as<NumericMatrix>(eqSS["ICMmat"]));
  Rcpp::NumericMatrix IDmat(Rcpp::as<NumericMatrix>(eqSS["IDmat"]));
  vector<double> age_brackets = Rcpp::as<vector<double> >(eqSS["age_brackets"]);
  vector<double> het_brackets = Rcpp::as<vector<double> >(eqSS["het_brackets"]);
  double ICM_Init = Rcpp::as<double>(eqSS["MaternalImmunity"]);
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // END: R -> C++ CONVERSIONS
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  Rcpp::Rcout << "Matrix unpacking working!\n";
  
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // START: HUMAN INITIALISATION FROM EQUILIBRIUM
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  // Use read in R equilibrium state to then allocate each individual accordingly given tehir age and biting heterogeneity
  
  // Initialisation temp variables required
  std::vector<int> infected_human_count;
  infected_human_count.reserve(parameters.g_N);
  int id_counter = 0;
  int num_age_brackets = age_brackets.size();
  int num_het_brackets = het_brackets.size();
  int age_bracket_in = 0;
  int het_bracket_in = 0;
  std::vector<double> infection_state_probability(6);
  
  // Temporary strain and pending strain vector needed for initialisation
  std::vector<Strain> temp_strains(100);
  std::vector<int> temp_infection_time_realisation_vector; temp_infection_time_realisation_vector.reserve(100);         
  std::vector<Person::InfectionStatus> temp_infection_state_realisation_vector; temp_infection_state_realisation_vector.reserve(100); 
  std::vector<boost::dynamic_bitset<> > temp_infection_barcode_realisation_vector; temp_infection_barcode_realisation_vector.reserve(100);   
  
  for (unsigned int n=0; n < parameters.g_N; n++) 
  {
    
    // Set their id
    population.push_back(Person(parameters));
    population[n].set_m_person_ID(id_counter++);
    
    // first find out what age bracket they are in
    // are they in the last - do this first for those people with ages that happen to be above the max age bracket
    if(population[n].get_m_person_age() >= (age_brackets[num_age_brackets-2]))
    {
      age_bracket_in = num_age_brackets-1;
    } 
    if(population[n].get_m_person_age() < age_brackets[0])
    {
      age_bracket_in = 0;
    }
    // if not then loop up to last-1
    for(int age_i = 1 ; age_i < (num_age_brackets-1) ; age_i++)
    {
      if(population[n].get_m_person_age() >= age_brackets[age_i-1] && population[n].get_m_person_age() < age_brackets[age_i])
      {
        age_bracket_in = age_i;
      }
    }
    
    // second find out what heterogeneity bracket they are in
    // are they in the last - do this first for those people with heterogeneities that happen to be above the max 
    if(population[n].get_m_individual_biting_rate() >= het_brackets[num_het_brackets-2]) 
    {
      het_bracket_in = num_het_brackets-1;
    }
    if(population[n].get_m_individual_biting_rate() < het_brackets[0])
    {
      het_bracket_in = 0;
    }
    for(int het_i = 1 ; het_i < (num_het_brackets-1) ; het_i++)
    {
      if(population[n].get_m_individual_biting_rate() >= het_brackets[het_i-1] && population[n].get_m_individual_biting_rate() < het_brackets[het_i])
      {
        het_bracket_in = het_i;
      }
    }
    
    infection_state_probability[0] = Smat(age_bracket_in, het_bracket_in);
    infection_state_probability[1] = Dmat(age_bracket_in, het_bracket_in);
    infection_state_probability[2] = Amat(age_bracket_in, het_bracket_in);
    infection_state_probability[3] = Umat(age_bracket_in, het_bracket_in);
    infection_state_probability[4] = Tmat(age_bracket_in, het_bracket_in);
    infection_state_probability[5] = Pmat(age_bracket_in, het_bracket_in);
    
    // Sample their infection state given the probabilities of being in any state given their age and het compartment
    population[n].set_m_infection_state(static_cast<Person::InfectionStatus>(sample1(infection_state_probability, std::accumulate(infection_state_probability.begin(), infection_state_probability.end(),0.0))));
    population[n].set_m_IB(IBmat(age_bracket_in, het_bracket_in));
    population[n].set_m_ICA(ICAmat(age_bracket_in, het_bracket_in));
    population[n].set_m_ICM_init(ICM_Init);
    if(population[n].get_m_person_age() == 0){
      population[n].set_m_ICM(ICM_Init); 
    }
    else
    {
    population[n].set_m_ICM(ICM_Init * exp(-population[n].get_m_person_age() / parameters.g_dCM));
    }
    population[n].set_m_ID(IDmat(age_bracket_in, het_bracket_in));
    
    // Schedule change for those who are not susceptible
    if (population[n].get_m_infection_state() != Person::SUSCEPTIBLE)
    {
      population[n].schedule_m_day_of_InfectionStatus_change(parameters);
    }
    
    // Check if mother and if so increase the maternal immunity sum and total number
    if (population[n].get_m_person_age() > 20 * 365 && population[n].get_m_person_age() < 21 * 365)
    {
      parameters.g_sum_maternal_immunity += population[n].get_m_ICA();
      parameters.g_total_mums++;
    }
    
    // If they are infected, i.e. not S or P, then assign their strains and next strain clearance date
    if (population[n].get_m_infection_state() != Person::SUSCEPTIBLE && population[n].get_m_infection_state() != Person::PROPHYLAXIS) 
    {
      // TODO: Think about how we can correctly initialise MOI for a given EIR. Presumably there is a rarefaction of MOI vs EIR, and the MOI is lognormal*age_dependency
      population[n].set_m_number_of_strains(static_cast<int>(population[n].get_m_individual_biting_rate()) + 1);
      population[n].set_m_number_of_realised_infections(population[n].get_m_number_of_strains());
      population[n].schedule_m_day_of_strain_clearance(parameters);
      
      // If they are treated the strains we allocate should have no state change
      if (population[n].get_m_infection_state() == Person::TREATED)
      {
        for (int s = 0; s < population[n].get_m_number_of_strains(); s++)
        {
          temp_strains.push_back( 
            Strain(
              Strain::generate_next_barcode(),
              Strain::m_transition_vector[static_cast<int>(population[n].get_m_infection_state())],
                                         0
            )
          );
        }
      }
      
      // If they are not then they will have state changes assumed equal to the human
      else
      {
        for (int s = 0; s < population[n].get_m_number_of_strains(); s++)
        {
          temp_strains.push_back(
            Strain(
              Strain::generate_next_barcode(),
              Strain::m_transition_vector[static_cast<int>(population[n].get_m_infection_state())],
                                         population[n].get_m_day_of_InfectionStatus_change()
            )
          );
        }
      }
      
      // Set infection time,state and barcode realisations
      for (int s = 0; s < population[n].get_m_number_of_strains(); s++)
      {
        temp_infection_time_realisation_vector.emplace_back(parameters.g_current_time);
        temp_infection_state_realisation_vector.emplace_back(population[n].get_m_infection_state());
        temp_infection_barcode_realisation_vector.emplace_back(temp_strains[s].get_m_barcode());
      }
      
      // Set human strains and associated times of strain state changing and clearing
      population[n].set_m_active_strains(temp_strains);
      population[n].set_m_day_of_next_strain_state_change();
      population[n].set_m_infection_time_realisation_vector_from_vector(temp_infection_time_realisation_vector);
      population[n].set_m_infection_state_realisation_vector(temp_infection_state_realisation_vector);
      population[n].set_m_infection_barcode_realisation_vector_from_vector(temp_infection_barcode_realisation_vector);
      temp_strains.clear();
      temp_infection_time_realisation_vector.clear();
      temp_infection_state_realisation_vector.clear();
      temp_infection_barcode_realisation_vector.clear();
      infected_human_count.push_back(n);
      
    }
    
    // Set the next event day
    population[n].set_m_day_of_next_event();
    
    // Add to the static zeta vector which is required for calculating the overall probability of being bitten, pi
    zeta_vector[n] = population[n].get_m_individual_biting_rate();
  }
  
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // END: HUMAN INITIALISATION FROM EQUILIBRIUM
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  Rcpp::Rcout << "Human initilisation working\n";
  
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // START: MOSQUITO INITIALISATION FROM EQUILIBRIUM
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  // Mosquito initialisation preallocation
  std::vector<double> mosquito_status_eq{ Sv, Ev, Iv };
  double mosquito_status_eq_sum = mosquito_status_eq[0] + mosquito_status_eq[1] + mosquito_status_eq[2];
  int parent_source = 0;
  std::vector<unsigned short int> pending_oocyst_time{ 0 };
  std::vector<boost::dynamic_bitset<> > pending_oocyst_barcode_male{ Strain::temp_barcode };
  std::vector<boost::dynamic_bitset<> > pending_oocyst_barcode_female{ Strain::temp_barcode };
  
  // exported barcodes vector set up
  if(parameters.g_spatial_type == Parameters::METAPOPULATION)
  {
    // parameters.g_spatial_export_counter = 0;
    // parameters.g_exported_barcodes.reserve(parameters.g_spatial_exports);
    // parameters.g_spatial_import_counter = 0;
    // parameters.g_imported_barcodes.reserve(parameters.g_spatial_imports);
  }
  
  
  Rcpp::Rcout << "Mosquito preallocation initilisation working\n";
  Rcpp::Rcout << scourge.capacity() << "\n";
  Rcpp::Rcout << infected_human_count.size() << "\n";
  
  // mosquito initialisation
  for (unsigned int n = 0; n < scourge.capacity(); n++)
  {
    // Set their id
    scourge.push_back(Mosquito(parameters));
    scourge[n].set_m_mosquito_ID(n);
    
    // Set if on/off season
    if (n >= parameters.g_scourge_today) scourge[n].set_m_mosquito_off_season(true);
    
    // Set infection status
    scourge[n].set_m_mosquito_infection_state(static_cast<Mosquito::InfectionStatus>(sample1(mosquito_status_eq, mosquito_status_eq_sum)));
    
    // Deal with if the mosquito is infected
    // Here we assume they have no pending oocysts but one already burst
    if (scourge[n].get_m_mosquito_infection_state() == Mosquito::INFECTED)
    {
      // Pick a random human for the source of the strain
      parent_source = runiform_int_1(0, infected_human_count.size() - 1);
      
      // If human has only one strain then assign the first strain
      if (population[infected_human_count[parent_source]].get_m_number_of_strains() == 1)
      {
        scourge[n].set_m_oocyst_barcode_male_vector(population[infected_human_count[parent_source]].get_m_person_strain_x(0).get_m_barcode());
        scourge[n].set_m_oocyst_barcode_female_vector(population[infected_human_count[parent_source]].get_m_person_strain_x(0).get_m_barcode());
      }
      // If human has more than one strain then pick random strain for both male and female mosquito barcode
      else
      {
        scourge[n].set_m_oocyst_barcode_male_vector(population[infected_human_count[parent_source]].get_m_person_strain_x(runiform_int_1(0, population[infected_human_count[parent_source]].get_m_number_of_strains() - 1)).get_m_barcode());
        scourge[n].set_m_oocyst_barcode_female_vector(population[infected_human_count[parent_source]].get_m_person_strain_x(runiform_int_1(0, population[infected_human_count[parent_source]].get_m_number_of_strains() - 1)).get_m_barcode());
      }

      // exported barcodes vector add barcode
      if(parameters.g_spatial_type == Parameters::METAPOPULATION)
      {
        // if(parameters.g_spatial_export_counter < parameters.g_spatial_exports)
        //   {
        //   parameters.g_exported_barcodes.emplace_back(population[infected_human_count[parent_source]].get_m_person_strain_x(0).get_m_barcode());
        //   parameters.g_spatial_export_counter++;
        // }
      }
      
      // Set the oocyst rupture count to 1
      scourge[n].set_m_ruptured_oocyst_count(1);
      scourge[n].m_mosquito_infected = true;

    }
    
    // Deal with if the mosquito is exposed
    // Here we will assume there is one pending oocyst with a random realisation time of up to delay.mos days for the purpose of initialisation
    if (scourge[n].get_m_mosquito_infection_state() == Mosquito::EXPOSED)
    {
      // Pick a random human for the source of the strain
      parent_source = runiform_int_1(0, infected_human_count.size() - 1);

      // If human has only one strain then assign the first strain
      if (population[infected_human_count[parent_source]].get_m_number_of_strains() == 1)
      {
        pending_oocyst_time[0] = static_cast<unsigned short int>(runiform_int_1(parameters.g_current_time, static_cast<int>(parameters.g_delay_mos)) + 1);
        pending_oocyst_barcode_male[0] = population[infected_human_count[parent_source]].get_m_person_strain_x(0).get_m_barcode();
        pending_oocyst_barcode_female[0] = population[infected_human_count[parent_source]].get_m_person_strain_x(0).get_m_barcode();
        scourge[n].set_m_oocyst_rupture_time_vector(pending_oocyst_time);
        scourge[n].set_m_oocyst_barcode_male_vector(pending_oocyst_barcode_male);
        scourge[n].set_m_oocyst_barcode_female_vector(pending_oocyst_barcode_female);

      }
      // If human has more than one strain then pick random strain for both male and female mosquito barcode
      else
      {

        pending_oocyst_time[0] = runiform_int_1(parameters.g_current_time, static_cast<int>(parameters.g_delay_mos)) + 1;
        pending_oocyst_barcode_male[0] = population[infected_human_count[parent_source]].get_m_person_strain_x(runiform_int_1(0, population[infected_human_count[parent_source]].get_m_number_of_strains() - 1)).get_m_barcode();
        pending_oocyst_barcode_female[0] = population[infected_human_count[parent_source]].get_m_person_strain_x(runiform_int_1(0, population[infected_human_count[parent_source]].get_m_number_of_strains() - 1)).get_m_barcode();
        scourge[n].set_m_oocyst_rupture_time_vector(pending_oocyst_time);
        scourge[n].set_m_oocyst_barcode_male_vector(pending_oocyst_barcode_male);
        scourge[n].set_m_oocyst_barcode_female_vector(pending_oocyst_barcode_female);
      }
      
    }

    // schedule mosquito's next day of event
    scourge[n].schedule_m_day_of_next_event();
    
  }
  
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // END: MOSQUITO INITIALISATION FROM EQUILIBRIUM
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  Rcpp::Rcout << "Mosquito initilisation working\n";
  

  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // BEGIN: LOGGING AND RETURN
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  
  // convert exproted barcodes to vector of vector of bool
  std::vector<SEXP>  Exported_Barcodes_Booleans(parameters.g_spatial_total_exported_barcodes);
  std::vector<SEXP>  Exported_Oocysts_Booleans(parameters.g_spatial_total_exported_oocysts);
  if(parameters.g_spatial_type == Parameters::METAPOPULATION)
  {
    
    // export barocdes
    for(unsigned int temp_status_iterator = 0; temp_status_iterator < parameters.g_spatial_total_exported_barcodes ; temp_status_iterator++)
    {
      Exported_Barcodes_Booleans[temp_status_iterator] = bitset_to_sexp(parameters.g_spatial_exported_barcodes[temp_status_iterator]);
    }
    
    // export oocysts
    for(unsigned int temp_status_iterator = 0; temp_status_iterator < parameters.g_spatial_total_exported_oocysts ; temp_status_iterator++)
    {
      Exported_Oocysts_Booleans[temp_status_iterator] = bitset_to_sexp(parameters.g_spatial_exported_oocysts[temp_status_iterator]);
    }
    
  }

  
  
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
    // Something like passing in a function name within the paramList which is the 
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
  Rcpp::Rcout << "S | D | A | U | T | P:\n" ;
  
  for (int element = 0; element < 6; element++) 
  {
    status_eq[element] /= (parameters.g_N);
    Rcpp::Rcout << status_eq[element] << " | ";
  }
  
  // Create Rcpp loggers list
  Rcpp::List Loggers = Rcpp::List::create(Rcpp::Named("S")=status_eq[0],Rcpp::Named("D")=status_eq[1],Rcpp::Named("A")=status_eq[2],
                                          Rcpp::Named("U")=status_eq[3],Rcpp::Named("T")=status_eq[4],Rcpp::Named("P")=status_eq[5],
                                          Rcpp::Named("InfectionStates")=Infection_States, Rcpp::Named("Ages")=Ages, 
                                          Rcpp::Named("IB")=IB,Rcpp::Named("ICA")=ICA,Rcpp::Named("ICM")=ICM,Rcpp::Named("ID")=ID);
  
  
  // Create universe ptr for memory-continuiation
  Rcpp::XPtr<Universe> universe_ptr(new Universe{ population, psi_vector, zeta_vector, pi_vector, scourge, parameters},
                                    true);
  
  // Return Named List with pointer and loggers
  // If spatial also required then export the barcodes
  if(parameters.g_spatial_type == Parameters::METAPOPULATION)
  {
  return Rcpp::List::create(Rcpp::Named("Ptr") = universe_ptr, Rcpp::Named("Loggers")=Loggers, 
                            Rcpp::Named("Exported_Barcodes")=Exported_Barcodes_Booleans,
                            Rcpp::Named("Exported_Oocysts")=Exported_Oocysts_Booleans);
  } 
  else
  {
  return Rcpp::List::create(Rcpp::Named("Ptr") = universe_ptr, Rcpp::Named("Loggers")=Loggers);
  }
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // fini
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
}


