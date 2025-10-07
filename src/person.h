//
//  magenta
//  person.h
//
//  Created: OJ on 12/01/2017
//  Most recent edits: OJ on 10/03/2017
//
//  Distributed under the MIT software licence - see Notes.c file for details
//
//  Class describing an individidual.
//
// ---------------------------------------------------------------------------


#ifndef PERSON_H
#define PERSON_H

#include "mosquito.h"
#include <vector>
#include <queue>
#include <deque>
#include <algorithm>
#include <numeric>
#include <cassert> // for error checking

using namespace std;

class Person {
  
public:
  enum InfectionStatus
  {
    SUSCEPTIBLE,  // 0
    DISEASED,   // 1
    ASYMPTOMATIC, // 2
    SUBPATENT,    // 3
    TREATED,    // 4
    PROPHYLAXIS,  // 5
    NUMBER_OF_STATES = 6
  };
  
  enum TreatmentOutcome
  {
    NOT_CLINICAL,  // 0
    SUCCESFULLY_TREATED,   // 1
    LPF, // 2
    NOT_TREATED, // 3
    NUMBER_OF_TO_STATES = 4
  };
  
  enum RecrudescenceOutcome
  {
    NOT_RECRUDESCING,  // 0
    RECRUDESCING_LPF,   // 1
    NUMBER_OF_RO_STATES = 2
  };
  
  // Infection transition options. 
  const static std::vector<InfectionStatus> m_transition_vector;
  
private:
  
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // Person identifiers, age, infection states, transition rates etc.
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  int m_person_ID;            // member variable ID - fine to be public as constant
  int m_person_age;           // Person's age
  unsigned int m_nmf_age_band = 0;      // what age bandare they in for nmf
  
  InfectionStatus m_infection_state;    // Infection Status enum
  InfectionStatus m_temp_infection_state;    // Infection Status temp enum
  
  TreatmentOutcome m_treatment_outcome = NOT_TREATED;
  RecrudescenceOutcome m_recrudescence_outcome = NOT_RECRUDESCING;
  bool m_slow_parasite_clearance_bool = false; // Flag for whether they are currently SPC
  int m_drug_choice = 0; // Flag for what drug they last received (i.e. so we know their prophylactic length)
  int m_drug_choice_time = 0; // Flag for what drug they last received (i.e. so we know their prophylactic length)
  
  // Person's age dependent biting rate (psi) - See Griffin 2010 S1 for this specific origin
  double m_age_dependent_biting_rate;
  // Person's symptom success rate (phi) - See Griffin 2010 S1 for this specific origin
  double m_symptom_success_rate;
  // Person's individual relative biting rate (zeta) due to mosquito biting heterogeneity following a lognormal distribution
  double m_individual_biting_rate;    // Set for life
  
  // Person's transition probabilities, i.e. the probabilty, given their acquired immunity and the treatment seeking rate, that they 
  // become diseased, treated or asymptomatic. Initialisation just for memory initialisation reasons, not used values
  std::vector<double> m_transition_probabilities = { 0.33, 0.33, 0.33 };
  
  // Sum of transition probabilities. Will always be 1 unless transitioning from DISEASED, as seems incorrect to transition to ASYMPTOMATIC due to an infection if already diseased
  double m_sum_transition_probabilities = 1;
  
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // Person's bite and strain counters and other temporary variables initialised for speed.
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  // Number of bites individual is receiving in this time step - only calculated for those being biten
  int m_number_of_bites = 0;
  
  // Number of succesful bites individual is receiving in this time step - only calculated for those being biten
  int m_number_of_succesful_bites = 0;
  
  // Number of strains an individual posseses - this can include multiple same strains
  int m_number_of_strains = 0;
  
  // Vector of strains
  std::vector<Strain> m_active_strains;
  std::vector<Strain> m_resistant_strains;
  std::vector<Strain> m_post_treatment_strains;
  
  // Temporary strain to be deleted
  int m_temp_strain_to_be_deleted;
  
  // Vector of associated active strain contribution to onwards infection and associate variables 
  std::vector<double> m_active_strain_contribution;
  std::vector<double> m_active_female_strain_contribution;
  int m_contribution_counter = 0;
  double m_contribution_sum = 0;
  double m_female_contribution_sum = 0;
  bool m_cA_counter = false;
  
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // IMMUNITY VARIABLES 
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  // Person's pre-erythrocytic immunity (IB), i.e. infection blocking
  double m_IB;                // infection-blocking immunity 
  double m_biting_success_rate = 0;     // Temporary biting success rate (b) dependent on age and pre-erythocytic immunity
  double m_IB_last_boost_time = runif1(0.0, 1.0);       // time that IB was last boosted
  int m_IB_last_calculated_time = 0;      // time IB was last calculated
  
  // Person's acquired disease blocking immunity 
  double m_ICM_init;
  double m_ICM;               // maternally acquired immunty (ICM)
  double m_ICA;               // exposure acquired immunity (ICA)
  double m_ICA_last_boost_time = runif1(0.0, 1.0);        // time that ICA was last boosted
  int m_I_C_D_CM_last_calculated_time = 0;    // time IC and ID was last calculated
  
  // Person's blood stage immunity (I.D), i.e. reducing probability of detectiona ond onwards contribution to infectiousness
  double m_ID;                // infection-blocking immunity 
  double m_ID_last_boost_time = runif1(0.0, 1.0);       // time that ID was last boosted
  double m_cA;                // individuals contribution to onwards infection from asymptomatic case
  
  // Important to note the above runifs as these are key in ensuring that the last boost times represent random times in the day
  // that biting occurs
  double m_immunity_boost_float = runif1(0.0, 1.0);
  
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // EVENT VARIABLES 
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  int m_day_of_InfectionStatus_change = 0;                                    // Person's time of InfectionStatus state change
  int m_day_of_strain_clearance = 0;                                          // Person's time of clearing a strain
  int m_day_of_death = 0;                                                     // Person's time of death
  int m_day_of_next_strain_state_change = std::numeric_limits<int>::max();    // Person's next strain state change day. Start with very large number due to person functions being written to compare against thsis often
  int m_day_of_next_event = 0;                                                // Person's closest event day
  bool m_more_than_one_strain_to_change_today_bool = false;                   // bool that declares whether there are more than one strain changing states today
  int m_day_last_treated = 0;                                                 // Day last treated
  int m_day_prophylaxis_wanes = 0;                                            // Day prophylaxis ends - this is just for tracking prophylaxis for recrudescent infections otherwise given by duration in P
  int m_temp_strain_to_next_change = 0;                                       // temp variable so we know what strain is changing next
  int m_day_of_nmf = std::numeric_limits<int>::max();                         // Day of nmf
  int m_temp_int = 0;                                                         // Temp integer for individual. This rather than a static int or new declaration so that it's thread safe and only one construct (?)
  
  std::vector<int> m_infection_time_realisation_vector{};         // First pending infection time in position 0 to handle multiple infections times that have not been realised yet
  std::vector<InfectionStatus> m_infection_state_realisation_vector{};  // First pending infection state in position 0 to handle multiple infections states that have not been realised yet
  std::vector<boost::dynamic_bitset<>> m_infection_barcode_realisation_vector{};    // First pending infection barcode in position 0 to handle multiple infections states that have not been realised yet
  
  int m_infection_realisation_empty_catch = 1;  // Variable that allows a check for empty vectors when dealing with more than one infection realisation on a day
  unsigned int m_number_of_realised_infections = 0;  // Count of realised infections
  int m_gametocytogenic_infections = 0;         // Count of any realised infections  occurred earlier than the current time - 12.5 days.
  std::vector<int> m_gametocytogenic_strains{}; // Vector of which strains are gametocytogenic - strains may not be always be in order of being acquired as we swap and pop when clearing strains 
  
public:
  
  // Default class constructor
  Person() {};
  
  // Non default class constructor
  Person(const Parameters &parameters);
  
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // GETTERS
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  // Get person's ID
  int get_m_person_ID() { return(m_person_ID); }
  
  // Get person's age
  int get_m_person_age() { return(m_person_age); }
  
  // Get person's infection status
  InfectionStatus get_m_infection_state() { return(m_infection_state); }
  
  // Get person's treatment outcome
  TreatmentOutcome get_m_treatment_outcome() { return(m_treatment_outcome); }

  // Get person's recrudescence outcome
  RecrudescenceOutcome get_m_recrudescence_outcome() { return(m_recrudescence_outcome); }

  // Get whether the person is currently experiencing slow parasite clearance
  bool get_m_slow_parasite_clearance_bool() { return(m_slow_parasite_clearance_bool); }

  // Get person's drug choice
  int get_m_drug_choice() { return(m_drug_choice); }

  // Get person's next non-malarial fever event day
  int get_m_day_of_nmf() { return(m_day_of_nmf); }

  // Get person's current non-malarial fever age band index
  unsigned int get_m_nmf_age_band() { return(m_nmf_age_band); }
  
  // Get person's age-dependent immunity
  double get_m_age_dependent_biting_rate() { return(m_age_dependent_biting_rate); }
  
  // Get person's number of strains
  int get_m_number_of_strains() { return(m_number_of_strains); }
  
  // Get person's temp strain to be next changed
  int get_m_temp_strain_to_next_change() { return(m_temp_strain_to_next_change); }
  
  // Get person's day of next strain state change
  int get_m_day_of_next_strain_state_change() { return(m_day_of_next_strain_state_change); }
  
  // Get person's number of realised infectons
  unsigned int get_m_number_of_realised_infections() { return(m_number_of_realised_infections); }
  
  // Get person's specific strain
  Strain get_m_person_strain_x(int x) { return(m_active_strains[x]); }
  
  // Get person's individual biting rate
  double get_m_individual_biting_rate() { return(m_individual_biting_rate); }
  
  // Get person's IB
  double get_m_IB() { return (m_IB); }
  
  // Get person's ICA
  double get_m_ICA() { return (m_ICA); }
  
  // Get person's ICM
  double get_m_ICM() { return (m_ICM); }
  
  // Get person's cA
  double get_m_cA() { return (m_cA); }
  
  // Get person's ICM_init
  double get_m_ICM_init() { return (m_ICM_init); }
  
  // Get person's ID
  double get_m_ID() { return (m_ID); }
  
  // Get person's IB_last_boost_time
  double get_m_IB_last_boost_time() { return (m_IB_last_boost_time); }
  
  // Get person's ICA_last_boost_time
  double get_m_ICA_last_boost_time() { return (m_ICA_last_boost_time); }
  
  // Get person's ID_last_boost_time
  double get_m_ID_last_boost_time() { return (m_ID_last_boost_time); }
  
  // Get person's IB_last_calculated_time
  int get_m_IB_last_calculated_time() { return (m_IB_last_calculated_time); }
  
  // Get person's I_C_D_CM_last_calculated_time
  int get_m_I_C_D_CM_last_calculated_time() { return (m_I_C_D_CM_last_calculated_time); }
  
  // Get person's immunity boost float
  double get_m_immunity_boost_float() { return (m_immunity_boost_float); }
  
  // Get person's time of InfectionStatus state change
  int get_m_day_of_InfectionStatus_change() { return (m_day_of_InfectionStatus_change); }
  
  // Get person's time of strain clearance
  int get_m_day_of_strain_clearance() { return (m_day_of_strain_clearance); }
  
  // Get person's time of death
  int get_m_day_of_death() { return (m_day_of_death); }
  
  // Get person's next event day
  int get_m_day_of_next_event() { return(m_day_of_next_event); }
  
  // Get person's day of being treated
  int get_m_day_last_treated() { return(m_day_last_treated); }
  
  // Get person's day of prophylaxis waning
  int get_m_day_prophylaxis_wanes() { return(m_day_prophylaxis_wanes); }
  
  // Get person's infection time realisation vector
  std::vector<int> get_m_infection_time_realisation_vector() { return (m_infection_time_realisation_vector); }

  // Get person's infection state realisation vector
  std::vector<InfectionStatus> get_m_infection_state_realisation_vector() { return (m_infection_state_realisation_vector); }
  
  // Get person's infection barcode realisation vector
  std::vector<boost::dynamic_bitset<>> get_m_infection_barcode_realisation_vector() { return (m_infection_barcode_realisation_vector); }
  
  // Get person's active strains vector
  std::vector<Strain> get_m_active_strains() { return (m_active_strains); }
  
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // SEMI-GETTERS
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  // Work out if reciprocal infection happened
  bool reciprocal_infection_boolean(const Parameters &parameters);
  
  // Work out if late parasitological failure happened
  bool late_paristological_failure_boolean(const Parameters &parameters);
  
  // Work out if reciprocal infection happened
  std::vector<boost::dynamic_bitset<>> sample_two_barcodes(const Parameters &parameters);
  
  // Work out chance of detecting malaria
  double q_fun(const Parameters &parms) {
    
    double fd =  1 - ((1-parms.g_fD0) / (1 + pow((m_person_age/parms.g_aD),parms.g_gD)));
    return(parms.g_d1 + ((1-parms.g_d1) / (1 + pow((m_ID/parms.g_ID0),parms.g_kD)*fd)));
    
  }
  
  // Malaria Detected by RDT/Micrscopy
  bool detectable_malaria(const Parameters &parms) {
    
    // are they high parasitaemia
    if (m_infection_state == DISEASED || m_infection_state == TREATED) {
      return(true);
    }
    
    // are they asymptomatic if so use q_fun
    if (m_infection_state == ASYMPTOMATIC) {
      return(rbernoulli1(q_fun(parms)));
    }
    
    return(false);
    
  }
  
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // SETTERS
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  // Set person ID
  void set_m_person_ID(int x) { m_person_ID = x; }
  
  // Get person's age
  void set_m_person_age(int x) { m_person_age = x; }
  
  // Set person's infection state
  void set_m_infection_state(InfectionStatus x) { m_infection_state = x; }
  
  // Set person's treatment outcome
  void set_m_treatment_outcome(TreatmentOutcome x) { m_treatment_outcome = x; }

  // Set person's recrudescence outcome
  void set_m_recrudescence_outcome(RecrudescenceOutcome x) { m_recrudescence_outcome = x; }

  // Set whether the person is currently experiencing slow parasite clearance
  void set_m_slow_parasite_clearance_bool(bool x) { m_slow_parasite_clearance_bool = x; }

  // Get person's drug choice
  void set_m_drug_choice(int x) { m_drug_choice = x; }

  // Set person's next non-malarial fever event day directly
  void set_m_day_of_nmf_from_saved(int x) { m_day_of_nmf = x; }

  // Set person's non-malarial fever age band index directly
  void set_m_nmf_age_band(unsigned int x) { m_nmf_age_band = x; }
  
  // Set person's individual biting rate
  void set_m_individual_biting_rate(double x) { m_individual_biting_rate = x; };
  
  // TODO: THESE ARE NOT FORMAL SETTERS AND SHOULD BE REWORDED TO SOMETHING LIKE CALCULATE
  // --------------------------------------------------------------------------------------------------------
  
  // Set person's individual biting rate
  double set_initial_m_individual_biting_rate(double zeta_meanlog, double zeta_sdlog);
  
  // Set person's relative chance of biting (psi)
  double set_m_age_dependent_biting_rate(double rho, double a0);
  
  // Set person's intial age given the maximum and average age
  int set_initial_m_person_age(double average_age, int max_age);
  
  // Set person's initial death day
  void set_initial_m_day_of_death(const Parameters &parameters);
  
  // Set person's transition probabilities
  void set_m_transition_probabilities(double pD, double pT, double pA) { m_transition_probabilities = { pD, pT, pA }; }
  
  // Set day of next strain state change
  void set_m_day_of_nmf(const Parameters &parameters);
  
  // Set day of next strain state change
  void set_m_day_of_next_strain_state_change();
  
  // Set day of next event
  void set_m_day_of_next_event();
  
  // --------------------------------------------------------------------------------------------------------
  
  // Set person's number of bites being received
  void set_m_number_of_bites(int x) { m_number_of_bites = x; }
  
  // Set person's number of bites being received
  void set_m_number_of_strains(int x) { m_number_of_strains = x; }
  
  // Set person's number of realised infectons
  void set_m_number_of_realised_infections(unsigned int x) { m_number_of_realised_infections = x; }
  
  // Set person's IB
  void set_m_IB(double x) { m_IB = x; }
  
  // Set person's ID
  void set_m_ID(double x) { m_ID = x; }
  
  // Set person's ICA
  void set_m_ICA(double x) { m_ICA = x; }
  
  // Set person's ICM
  void set_m_ICM(double x) { m_ICM = x; }
  
  // Set person's ICM_init
  void set_m_ICM_init(double x) { m_ICM_init = x; }
  
  // Set person's IB_last_boost_time
  void set_m_IB_last_boost_time(double x) { m_IB_last_boost_time = x; }
  
  // Set person's ICA_last_boost_time
  void set_m_ICA_last_boost_time(double x) { m_ICA_last_boost_time = x; }
  
  // Set person's ID_last_boost_time
  void set_m_ID_last_boost_time(double x) { m_ID_last_boost_time = x; }
  
  // Set person's IB_last_calculated_time
  void set_m_IB_last_calculated_time(int x) { m_IB_last_calculated_time = x; }
  
  // Set person's I_C_D_CM_last_calculated_time
  void set_m_I_C_D_CM_last_calculated_time(int x) { m_I_C_D_CM_last_calculated_time = x; }
  
  // Set person's immunity boost float
  void set_m_immunity_boost_float(double x) { m_immunity_boost_float = x; }
  
  // Set person's time of InfectionStatus state change
  void set_m_day_of_InfectionStatus_change(int x) { m_day_of_InfectionStatus_change = x; }
  
  // Set person's time of strain clearance
  void set_m_day_of_strain_clearance(int x) { m_day_of_strain_clearance = x; }
  
  // Set person's time of death
  void set_m_day_of_death(int x) { m_day_of_death = x; }
  
  // Set person's day of being treated
  void set_m_day_last_treated(int x) { m_day_last_treated = x; }
  
  // Set person's next event day
  void set_m_day_prophylaxis_wanes(int x) { m_day_prophylaxis_wanes = x; }
  
  // Set person's infection time realisation vector
  void set_m_infection_time_realisation_vector_from_vector(std::vector<int> x) 
  {
    for (unsigned int i = 0; i < x.size(); i++) 
    {
      m_infection_time_realisation_vector.emplace_back(x[i]);
    }
  }
  
  // Set person's infection state realisation vector
  void set_m_infection_state_realisation_vector_from_vector(std::vector<int> x)
  {
    for (unsigned int i = 0; i < x.size(); i++) {
      m_infection_state_realisation_vector.emplace_back(static_cast<Person::InfectionStatus>(x[i]));
    }
  }
  
  // Set person's infection state realisation vector
  void set_m_infection_state_realisation_vector(std::vector<InfectionStatus> x) { m_infection_state_realisation_vector = x; }
  
  // Set person's infection barcode realisation vector
  void set_m_infection_barcode_realisation_vector_from_vector(std::vector<boost::dynamic_bitset<>> x)
  {
    for (unsigned int i = 0; i < x.size(); i++) {
      m_infection_barcode_realisation_vector.emplace_back(x[i]);
    }
  }
  
  // Set person's barcode realisation vector from a vector<vector<bool> >
  void set_m_infection_barcode_realisation_vector_from_vector_of_vector_bool(std::vector<std::vector<bool> > x)
  {
    
    m_infection_barcode_realisation_vector.reserve(x.size());
    unsigned int temp_barcode_iterator = 0;
    
    for (unsigned int i = 0; i < x.size(); i++) {
      
      // fetch vector<bool> and turn into barcode
      for (temp_barcode_iterator = 0; temp_barcode_iterator < Parameters::g_barcode_length; temp_barcode_iterator++)
      {
        Strain::temp_barcode[temp_barcode_iterator] = x[i][temp_barcode_iterator];
      }
      m_infection_barcode_realisation_vector.emplace_back(Strain::temp_barcode);
    }
  }
  
  // Set person's active strains vector
  void set_m_active_strains(std::vector<Strain> x) { m_active_strains = x; }
  
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // ALLOCATIONS
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  // Allocate bite to person
  void allocate_bite(Parameters &parameters, Mosquito &mosquito);
  
  // Allocate an infection to person
  void allocate_infection(Parameters &parameters, Mosquito &mosquito);
  
  // Allocate a strain to a person
  void allocate_strain_with_push(Strain x) { m_active_strains.push_back(x); }
  
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // SCHEDULERS & DRAWS
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  // Draw day of infection state change 
  int draw_m_day_of_InfectionStatus_change(const Parameters &parameters);
  
  // Schedule person's infection state change
  void schedule_m_day_of_InfectionStatus_change(const Parameters &parameters);
  
  // Schedule person's death day
  void schedule_m_day_of_death(const Parameters &parameters);
  
  // Schedule person's next strain clearance
  // DEPRECATED - NOT USED ANY MORE
  void schedule_m_day_of_strain_clearance(const Parameters &parameters);
  
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // UPDATERS - functions that happen at the end of a day, i.e. ageing, dying, clearing strains etc. 
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  // Clear one strain from an individual
  void individual_strain_clearance();
  
  // Clear all strains from an individual
  void all_strain_clearance();
  
  // Clear the last strain if it would have been cleared by prophylaxis
  void clear_strain_if_prophylactic(const Parameters &parameters);
  
  // Recover to being susceptible, i.e. clearing all infections and strains and associated timings
  void recover(const Parameters &parameters);
  
  // Treatment outcomes
  void treatment_outcome(const Parameters &parameters);
  
  // Late parasitological failure
  void late_paristological_failure(const Parameters &parameters); 
  
  // Slow parasite clearance
  void slow_treatment_clearance(const Parameters &parameters);
  
  // Seek treatment for nmf
  void seek_nmf_treatment(const Parameters &parameters);
  
  // Kill person, i.e. reset age to 0, infections to 0, state to susceptible, immunities reset etc
  void die(const Parameters &parameters);
  
  // Daily update function, i.e. increase ages, set maternal booleans, calculate biting proportions etc
  double update(Parameters &parameters);
  
  // Event handle, i.e. if any of person's death, state change, strain clearance or state change realisation days are today, respond accordingly
  void event_handle(const Parameters &parameters);
  
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // LOGGERS - functions that report summaries on an individual
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  // Caluclate daily incidence for population and 0-5 years, i.e. would they cause an incident case today.
  // Returns 2 if 0-5 and incident, 1 if just incident, and 0 if neither
  int log_daily_incidence(const Parameters &parameters);
  
  // Update immunities so they are correct for today - helpful when logging to know actual immunities
  void update_immunities_to_today(const Parameters &parameters);
  
};

#endif // PERSON_H