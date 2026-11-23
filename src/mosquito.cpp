#include "mosquito.h"
#include "person.h"

// Non class constructor which will inialise a random strain
// First blood meal is the days until next feed when feeding each day with probability 1/(mean biting interval)
Mosquito::Mosquito(const Parameters &parameters) :
  m_day_of_death(rexpint1(parameters.g_mean_mosquito_age) + parameters.g_current_time ),
  m_day_of_next_blood_meal(static_cast<int>(std::log(unif_rand()) / std::log(1.0 - 1.0 / parameters.g_mosquito_next_biting_day_vector[0])) + parameters.g_current_time)
{
  // Reserve vector memory
  m_oocyst_rupture_time_vector.reserve(10);
  m_oocyst_barcode_male_vector.reserve(10);
  m_oocyst_barcode_female_vector.reserve(10);
  m_generated_sporozoites_vector.reserve(10);
  m_oocyst_remaining_spz_count.reserve(10);
  
}

// temporary mosquito variables
std::vector<boost::dynamic_bitset<> > gam_sampled(2);
int m_oocyst_pick = 0;

// ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// CHECKERS
// ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

// Function to quickly check if mosquito could carry resistant sporozoites
bool Mosquito::check_resistance(){
  
  // check the spz already we know about
  for(auto b : m_generated_sporozoites_vector){
    if(b.any()) {
      return(true);
    }
  }
  
  // check the females
  for(unsigned int f(0) ; f < m_ruptured_oocyst_count ; f++){
    if(m_oocyst_barcode_female_vector[f].any()) {
      return(true);
    }
  }
  
  // check the males
  for(unsigned int m(0) ; m < m_ruptured_oocyst_count ; m++){
    if(m_oocyst_barcode_male_vector[m].any()) {
      return(true);
    }
  }
  
  return(false);
}

// ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// ALLOCATIONS
// ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

// Allocate male and female gameteocyte resulting froom biting infected human
void Mosquito::allocate_gametocytes(const Parameters &parameters,
                                    boost::dynamic_bitset<> &gam_1,
                                    boost::dynamic_bitset<> &gam_2)
{
  // Push gametocyte barcodes and time of bursting and spz count remaining, i.e. 4
  m_oocyst_barcode_female_vector.emplace_back(gam_1);
  m_oocyst_barcode_male_vector.emplace_back(gam_2);
  // m_oocyst_barcode_female_vector.push_back(gam_1);
  // m_oocyst_barcode_male_vector.push_back(gam_2);
  m_oocyst_rupture_time_vector.emplace_back(parameters.g_current_time + parameters.g_mosquito_current_eip);
  m_oocyst_remaining_spz_count.emplace_back(4);
  
  // Update infection status if mosquito is susceptible
  if (m_mosquito_infection_state == SUSCEPTIBLE)
  {
    m_mosquito_infection_state = EXPOSED;
  }
  
}

// Handle bite, i.e. allocating gametocytes and other housekeeping
void Mosquito::handle_bite(Parameters &parameters, Person &person)
{
  
  parameters.g_total_mosquito_infections++;
  
  // EIP for all oocysts from this bite, drawn once per bite so pending oocysts stay in rupture order
  parameters.g_mosquito_current_eip = parameters.g_mosquito_eip_vector[parameters.g_mosquito_eip_counter++];
  if (parameters.g_mosquito_eip_counter == parameters.g_max_mosquito_biting_counter) parameters.g_mosquito_eip_counter = 0;
  
  // are we doing vector adaptation work. If so alter the oocyst numbers 
  if (parameters.g_gametocyte_sterilisation_flag){
    if (person.get_m_infection_state() == Person::TREATED || 
        (person.get_m_infection_state() == Person::ASYMPTOMATIC && person.get_m_treatment_outcome() == Person::LPF && person.get_m_day_last_treated() > parameters.g_current_time - 10) ) {
      parameters.g_oocyst_frequencies[parameters.g_oocyst_frequencies_counter] = std::ceil(parameters.g_oocyst_frequencies[parameters.g_oocyst_frequencies_counter]*parameters.g_oocyst_reduction_by_artemisinin);
    }
  }
  
  // if we are doing spatial then use the imported oocysts first 
  if(parameters.g_spatial_imported_mosquito_infection_counter < parameters.g_spatial_total_imported_mosquito_infections ||
     parameters.g_percentage_imported_mosquito_infections == 1.0)
  {
    
    // handle if metapopulation
    if(parameters.g_spatial_type == Parameters::METAPOPULATION)
    {
      // loop ove the multiple oocyst frequencies
      for(int o_i = parameters.g_spatial_imported_oocyst_frequencies[parameters.g_spatial_imported_mosquito_infection_counter] ;
          o_i < parameters.g_spatial_imported_oocyst_frequencies[parameters.g_spatial_imported_mosquito_infection_counter+1] ;
          o_i ++)
      {
        allocate_gametocytes(parameters, parameters.g_spatial_imported_oocysts[o_i*2],parameters.g_spatial_imported_oocysts[(o_i*2)+1]);
      }
      // also need to export oocysts here
      // TODO in finishing metapopulation
    } 
    // else for island simulations
    else 
    {
      // allocate imported gametocytes
      for(int o_i = 0; o_i < parameters.g_oocyst_frequencies[parameters.g_oocyst_frequencies_counter]; o_i++)
      {  
        gam_sampled[0] = Strain::generate_next_barcode();
        gam_sampled[1] = Strain::generate_next_barcode();
        
        // are we doing vector adaptation
        if (parameters.g_vector_adaptation_flag){
          
          // check if the oocyst would have been produced by judging the female gametocyte for adaptation
          
          // Check if all vector adapted loci are true
          if (Strain::all_at_positions(gam_sampled[0], parameters.g_vector_adaptation_loci)) {
            allocate_gametocytes(parameters, gam_sampled[0], gam_sampled[1]);
          } else {
            // if not would it have still made it through
            if(rbernoulli1(parameters.g_local_oocyst_advantage)) {
              allocate_gametocytes(parameters, gam_sampled[0], gam_sampled[1]);
            }
          }
        }
        else 
        {
          allocate_gametocytes(parameters, gam_sampled[0], gam_sampled[1]);
        }
      }
    }
    
    // and increase the counter
    parameters.g_spatial_imported_mosquito_infection_counter++;
    
  }
  else 
  {
    // // multiple oocyst handling
    for(int oocyst_freq = 0 ; oocyst_freq < parameters.g_oocyst_frequencies[parameters.g_oocyst_frequencies_counter] ; oocyst_freq++)
    {  
      gam_sampled = person.sample_two_barcodes(parameters);
      
      // are we doing vector adaptation
      if (parameters.g_vector_adaptation_flag){
        
        // check if the oocyst would have been produced by judging the female gametocyte for adaptation
        
        // Check if all vector adapted loci are true
        if (Strain::all_at_positions(gam_sampled[0], parameters.g_vector_adaptation_loci)) {
          allocate_gametocytes(parameters, gam_sampled[0], gam_sampled[1]);
        } else {
          // if not would it have still made it through
          if(rbernoulli1(parameters.g_local_oocyst_advantage)) {
            allocate_gametocytes(parameters, gam_sampled[0], gam_sampled[1]);
          }
        }
      }
      else 
      {
        allocate_gametocytes(parameters, gam_sampled[0], gam_sampled[1]);
      }
    }
  }
  
  // Increase oocyst counter and catch for overflow
  if(++parameters.g_oocyst_frequencies_counter == parameters.g_oocyst_frequencies_size) parameters.g_oocyst_frequencies_counter = 0;
  
  // Schedule next event
  schedule_m_day_of_next_event();
  
}

// Sample one sporozoite to be passed on to the human
boost::dynamic_bitset<> Mosquito::sample_sporozoite() {
  
  // are there any spz in the mosquito
  if(m_generated_sporozoite_count){
    // are we using a spz
    if(rbernoulli1(m_generated_sporozoite_count / (4*m_ruptured_oocyst_count) )) 
    {
      return(m_generated_sporozoites_vector[runiform_int_1(1, m_generated_sporozoite_count - 1)]);
    }
  }
  
  // sample from avaialble oocysts and then store the spz
  // does the mosquito have more than one oocyst to pick from 
  if (m_ruptured_oocyst_count == 1)
  {
    Strain::temp_barcode_pair = Strain::generate_recombinant_barcode(
      get_m_oocyst_barcode_male_vector(0),
      get_m_oocyst_barcode_female_vector(0)
    );
    
    // decrease the remaining  spz from this oocyst
    m_oocyst_remaining_spz_count[0] = m_oocyst_remaining_spz_count[0] - 2;
  } 
  else 
  {
    m_oocyst_pick = sample1_ints(m_oocyst_remaining_spz_count, (4*m_ruptured_oocyst_count)-m_generated_sporozoite_count);
    Strain::temp_barcode_pair = Strain::generate_recombinant_barcode(
      m_oocyst_barcode_male_vector[m_oocyst_pick],
                                  m_oocyst_barcode_female_vector[m_oocyst_pick]
    );
    
    // decrease the number of remaining spz from this oocyst
    m_oocyst_remaining_spz_count[m_oocyst_pick] = m_oocyst_remaining_spz_count[m_oocyst_pick] - 2;
  }
  
  m_generated_sporozoite_count = m_generated_sporozoite_count + 2;
  m_generated_sporozoites_vector.emplace_back(Strain::temp_barcode_pair[0]);
  m_generated_sporozoites_vector.emplace_back(Strain::temp_barcode_pair[1]);
  
  // which one was then passed on
  if(rbernoulli1(0.5)) {
    return(Strain::temp_barcode_pair[0]);
  } else {
    return(Strain::temp_barcode_pair[1]);
  }
  
  return(Strain::temp_barcode);
  
} 

// Allocate sporozoite arising from importation
void Mosquito::allocate_imported_sporozoite(boost::dynamic_bitset<> x) {
  m_generated_sporozoite_count++;
  m_generated_sporozoites_vector.emplace_back(x);
}

// ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// SCHEDULERS 
// ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


// Schedule mosquito's death day
// Also return true if the scheduled day is one such that we care about bites this mosquito takes, i.e.
// it will live long enought to become infectious
bool Mosquito::schedule_m_day_of_death(const Parameters &parameters)
{
  
  // Exponential waiting time plus current day and 1 so not the same day
  // If less than 11 also move their blood meal day  to after their death
  // and set next event to their death day as 
  // and return true so we don't bother resetting anything else about the mosquito in the die function
  m_day_of_death = rexpint1(parameters.g_mean_mosquito_age);
  
  // if (m_day_of_death < 11) {
  //   m_day_of_next_blood_meal = parameters.g_current_time + 12;
  //   
  //   m_day_of_death += parameters.g_current_time + 1;
  //   m_day_of_next_event = m_day_of_death;
  //   return(false);
  // }
  
  m_day_of_death += parameters.g_current_time + 1;
  return(true);
  
}

// Schedule mosquito's day of next blood meal
void Mosquito::schedule_m_day_of_blood_meal(Parameters &parameters)
{
  
  // Current day plus next biting. Next biting is a vector of values giving the next biting day. 
  // If there are no interventions this will always be 3, but will sometimes be 4 given interventions
  m_day_of_next_blood_meal = parameters.g_current_time + parameters.g_mosquito_next_biting_day_vector[parameters.g_mosquito_biting_counter++];
  
  // If we are at end of biting day vector with our counter then reset it to 0
  if(parameters.g_mosquito_biting_counter == parameters.g_max_mosquito_biting_counter) parameters.g_mosquito_biting_counter=0;
}

// Set day of next event
void Mosquito::schedule_m_day_of_next_event()
{
  
  // Start next event as death
  m_day_of_next_event = m_day_of_death;
  
  // Catch if queue is empty first
  // Check if the next pending oocyst rupture is the next event
  if (m_oocyst_rupture_time_vector.size() > m_ruptured_oocyst_count)
  {
    if (m_day_of_next_event > m_oocyst_rupture_time_vector[m_ruptured_oocyst_count] && m_oocyst_rupture_time_vector[m_ruptured_oocyst_count] > 0) {
      m_day_of_next_event = m_oocyst_rupture_time_vector[m_ruptured_oocyst_count];
    }
  }
  
  //}
  
  // Compare against strain clearance day
  if (m_day_of_next_event > m_day_of_next_blood_meal && m_day_of_next_blood_meal > 0) {
    m_day_of_next_event = m_day_of_next_blood_meal;
  }
  
}

// ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// UPDATERS 
// ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

// Kill mosquito, i.e. reset age to 0, infections to 0, state to susceptible, ino barcodes etc
// Returns true if mosquito lives long enough to surive incubation, otherwise false
bool Mosquito::die(Parameters & parameters)
{
  
  // Reset whether mosquito is off season to false, as if we are killing th mosquito we must be considering it
  m_mosquito_off_season = false;
  
  // Schedule death day. If returns true, i.e. mosquito won't survive incubation, then return early as
  // no need to reset everything else
  if (!schedule_m_day_of_death(parameters)) {
    return (false);
  }
  
  // Reset mosquito's age and infection status changes
  m_ruptured_oocyst_count = 0;
  m_generated_sporozoite_count = 0;
  m_day_of_next_event = 0;
  m_mosquito_infected = false;
  
  // Clear vectors
  m_oocyst_rupture_time_vector.clear();
  m_oocyst_barcode_male_vector.clear();
  m_oocyst_barcode_female_vector.clear();
  m_generated_sporozoites_vector.clear();
  m_oocyst_remaining_spz_count.clear();
  
  // Make mosquito susceptible again
  m_mosquito_infection_state = SUSCEPTIBLE;
  
  if (m_day_of_next_blood_meal == parameters.g_current_time)
  {
    
    // Schedule next blood meal
    schedule_m_day_of_blood_meal(parameters);
    
    // Set to biting today
    m_mosquito_biting_today = true;
    
  }
  
  // Schedule next event
  schedule_m_day_of_next_event();
  
  // Schedule next blood meal
  // schedule_m_day_of_blood_meal(parameters);
  // schedule_m_day_of_blood_meal(parameters);
  
  return(m_mosquito_biting_today);
}


// Daily update function, i.e. check for sesonal mosquitoes and then handle any events
bool Mosquito::update(Parameters &parameters)
{
  
  // assume not biting
  m_mosquito_biting_today = false;
  
  // if the mosquito is off season, and there is a positive deficit, i.e. a surplus of mosquitoes then return early
  if (m_mosquito_off_season && parameters.g_mosquito_deficit > 0) {
    return(false);
  } 
  
  // if the mosquito is off season, and there is a deficit, i.e. a need of mosquitoes then return true early as will be biting today
  // DEPRECATED IF STATEMENT
  // This will never be reached as when a new mosquito first comes back on season we don't update it we simply copy an already on season one
  // in order to keep the correct mosquito popualtion infection. 
  if (m_mosquito_off_season && parameters.g_mosquito_deficit < 0) {
    parameters.g_mosquito_deficit++;
    return(die(parameters));
  }
  
  // Handle any events that are happening today
  if (m_day_of_next_event == parameters.g_current_time) {
    // if event handle returns true the mosquito is biting today 
    return(event_handle(parameters)); 
  }
  
  // if there was no event or seasonality then return false
  return(false);
  
}

// Event handle, i.e. mosquito death, state change, biting handle.
// If it returns true then the mosquito is biting today
bool Mosquito::event_handle(Parameters &parameters)
{
  
  // If death occurs then no need to explore the other events
  if (m_day_of_death == m_day_of_next_event)
  {
    return(die(parameters));
  }
  else
  {
    // All other events could happen theoretically on the same day though so within same else block
    // Update blood meal day
    if (m_day_of_next_blood_meal == m_day_of_next_event)
    {
      
      // Schedule next blood meal
      schedule_m_day_of_blood_meal(parameters);
      
      // Set to biting today
      m_mosquito_biting_today = true;
      
    }
    
    // Realise oocyst rupturing if time to do so
    // If there are more pending oocyst times than ruptures then see if the pending is happening today
    if (m_oocyst_rupture_time_vector.size() > m_ruptured_oocyst_count)
    {
      
      // if the vector is not for today, i.e. oocyst does not burst today
      if (m_oocyst_rupture_time_vector[m_ruptured_oocyst_count] == m_day_of_next_event)
      {
        
        // Update the catch - this is needed for when two oocysts might burst on the same day
        m_oocyst_realisation_empty_catch = 1;
        
        // Make infection status infected
        m_mosquito_infection_state = INFECTED;
        m_mosquito_infected = true;
        
        while (m_oocyst_realisation_empty_catch == 1)
        {
          
          // Update the burst oocyst count
          m_ruptured_oocyst_count++;
          
          // Update the catching - i.e. if there are no oocysts to rupture today then do not update the catch
          if(m_oocyst_rupture_time_vector.size() > m_ruptured_oocyst_count)
          {
            if (m_oocyst_rupture_time_vector[m_ruptured_oocyst_count] != parameters.g_current_time)
            {
              m_oocyst_realisation_empty_catch = 0;
            }
          }
          else 
          {
            m_oocyst_realisation_empty_catch = 0;
          }
        }
        
        
      }
    }
    
    // Update the day of next event
    schedule_m_day_of_next_event();
    
    // Throw if the next event is in the past unless mosquito is off season
    assert((m_day_of_next_event > parameters.g_current_time || m_mosquito_off_season) &&
      "Mosquitoes stuck in the past");
    
    // If mosquito is biting today return true
    return(m_mosquito_biting_today);
  }
  
  
  
}
