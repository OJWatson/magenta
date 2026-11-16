#include "universe_utils.h"


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

// q function for detectability
double q_fun(const Parameters &parms, unsigned int age, double ID) {
  
  double fd =  1 - ((1-parms.g_fD0) / (1 + pow((age/parms.g_aD),parms.g_gD)));
  return(parms.g_d1 + ((1-parms.g_d1) / (1 + pow((ID/parms.g_ID0),parms.g_kD)*fd)));
}


//' Returns the population's parasite genetics for ibd style summarised by pibd for given sample size and state
//'
//' @param param_list param_list containing statePtr, sample_size, and sample_states
//' @return list of population information
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List population_get_genetics_ibd_df_n(Rcpp::List param_list) {
  
  // create universe and do parameter conversions
  Rcpp::XPtr<Universe> u_ptr = Rcpp::as<Rcpp::XPtr<Universe> > (param_list["statePtr"]);
  unsigned int sample_size = Rcpp::as<unsigned int>(param_list["sample_size"]);
  std::vector<unsigned int> sample_states = Rcpp::as<std::vector<unsigned int> >(param_list["sample_states"]);
  
  unsigned int pop_size = u_ptr->population.size();
  // handle if we are subsampling
  std::vector<int> pop_sample(pop_size);  
  std::iota(pop_sample.begin(), pop_sample.end(), 0);
  std::vector<unsigned int> all_states = {0, 1, 2, 3, 4, 5};
  
  // if we are subsampling states
  if(all_states != sample_states){
    pop_sample.clear();
    for(auto p : u_ptr->population) {
      if(std::find(sample_states.begin(), sample_states.end(), static_cast<unsigned int>(p.get_m_infection_state())) != sample_states.end()) {
        pop_sample.emplace_back(p.get_m_person_ID());
      }
    }
  }
  
  // if we are choosing a sample
  if(sample_size) {
    pop_sample = sample_without_replacement(pop_sample, sample_size);
  }
  pop_size = pop_sample.size();
  
  
  std::vector<double> pibd_d(pop_size, 0);
  std::vector<double> pibd_within_d(pop_size, 0);
  std::vector<double> pibd(pop_size, 0);
  std::vector<double> pibd_within(pop_size, 0);
  std::vector<unsigned int> ages(pop_size, 0);
  std::vector<unsigned int> states(pop_size, 0);
  std::vector<std::vector<unsigned int> > barcode_states(pop_size);
  std::vector<std::vector<boost::dynamic_bitset<> > > bitsets(pop_size);
  std::vector<std::vector<boost::dynamic_bitset<> > > bitsets_d(pop_size);
  

  // temps
  double q_i;
  unsigned int total_barcodes = 0;
  unsigned int el_i;
  std::vector<Strain> strains_i;
  std::vector<boost::dynamic_bitset<> > bitsets_i;
  std::vector<boost::dynamic_bitset<> > bitsets_d_i;
  std::vector<unsigned int> barcode_states_i;
  
  for (unsigned int el = 0; el < pop_size; el++) {
    
    // who are we sampling
    el_i = static_cast<unsigned int>(pop_sample[el]);
    
    // simple assignments first
    ages[el] = u_ptr->population[el_i].get_m_person_age();
    states[el] = static_cast<unsigned int>(u_ptr->population[el_i].get_m_infection_state());
    
    // what's their disease detectability
    q_i = q_fun(u_ptr->parameters, ages[el], u_ptr->population[el_i].get_m_ID());
    
    // clear temps
    bitsets_i.clear();
    bitsets_d_i.clear();
    strains_i.clear();
    barcode_states_i.clear();
    
    // get their strains
    strains_i = u_ptr->population[el_i].get_m_active_strains();
    
    // reserve space
    bitsets_i.reserve(strains_i.size());
    barcode_states_i.reserve(strains_i.size());
    barcode_states[el].reserve(strains_i.size());
    bitsets[el].reserve(strains_i.size());

    
    // and then their barcode bitsets for all D and those detected in A
    for(auto s : strains_i) {
      if(s.get_m_strain_infection_status() == Strain::DISEASED || s.get_m_strain_infection_status() == Strain::TREATED) {
        bitsets_i.emplace_back(s.get_m_barcode());
        bitsets_d_i.emplace_back(s.get_m_barcode());
        barcode_states_i.emplace_back(static_cast<unsigned int>(s.get_m_strain_infection_status()));
      } 
      if(s.get_m_strain_infection_status() == Strain::ASYMPTOMATIC){
        if(rbernoulli1(q_i)) {
          bitsets_i.emplace_back(s.get_m_barcode());
          barcode_states_i.emplace_back(static_cast<unsigned int>(s.get_m_strain_infection_status()));
        }  
      }
      total_barcodes++;
    }
    
    // store these for later
    bitsets[el] = bitsets_i;
    bitsets_d[el] = bitsets_d_i;
    barcode_states[el] = barcode_states_i;
    
    // whats the pibd within the individuals bitsets that could be detected and within just those at high parasitaemia
    pibd_within[el] = Strain::ibd_distance_mean_within_bitsets(bitsets_i, true);
    pibd_within_d[el] = Strain::ibd_distance_mean_within_bitsets(bitsets_d_i, true);
    
  }
  
  // now work out the pibd within the population
  double pibd_i = 0.0;
  double pibd_d_i = 0.0;
  std::vector<int> contacts(pop_size,0);
  std::vector<int> contacts_d(pop_size,0);
  std::vector<boost::dynamic_bitset<> > temp_bitset_vec;
  std::vector<boost::dynamic_bitset<> > temp_bitset_vec_d;
  std::vector<boost::dynamic_bitset<> > temp_bitset_vec2;
  std::vector<boost::dynamic_bitset<> > temp_bitset_vec_d2;
  
  
  for(unsigned int i = 0; i < (pop_size-1); i++) {
    if(bitsets[i].size()){
          temp_bitset_vec = bitsets[i];
          std::sort( temp_bitset_vec.begin(), temp_bitset_vec.end() );
          temp_bitset_vec.erase( unique( temp_bitset_vec.begin(), temp_bitset_vec.end() ), temp_bitset_vec.end() );
          
          temp_bitset_vec_d = bitsets_d[i];
          std::sort( temp_bitset_vec_d.begin(), temp_bitset_vec_d.end() );
          temp_bitset_vec_d.erase( unique( temp_bitset_vec_d.begin(), temp_bitset_vec_d.end() ), temp_bitset_vec_d.end() );
      for(unsigned int j = i+1; j < pop_size; j++) {
        if(bitsets[j].size()) {
          pibd_i = pibd_d_i = 0;
          
          temp_bitset_vec2 = bitsets[j];
          std::sort( temp_bitset_vec2.begin(), temp_bitset_vec2.end() );
          temp_bitset_vec2.erase( unique( temp_bitset_vec2.begin(), temp_bitset_vec2.end() ), temp_bitset_vec2.end() );
          for(unsigned int k1 = 0; k1 < temp_bitset_vec.size(); k1++){
            pibd_i = pibd_i + Strain::ibd_distance_of_bitset_a_and_x(temp_bitset_vec[k1], temp_bitset_vec2.begin(), temp_bitset_vec2.end());
          }
          
          temp_bitset_vec_d2 = bitsets_d[j];
          std::sort( temp_bitset_vec_d2.begin(), temp_bitset_vec_d2.end() );
          temp_bitset_vec_d2.erase( unique( temp_bitset_vec_d2.begin(), temp_bitset_vec_d2.end() ), temp_bitset_vec_d2.end() );
          for(unsigned int k2 = 0; k2 < temp_bitset_vec_d.size(); k2++){
            pibd_d_i = pibd_d_i + Strain::ibd_distance_of_bitset_a_and_x(temp_bitset_vec_d[k2], temp_bitset_vec_d2.begin(), temp_bitset_vec_d2.end());
          }
          
          // running total of the number of bitset comparisons we have made for each individual
          contacts[i] += temp_bitset_vec.size() * temp_bitset_vec2.size();
          contacts[j] += temp_bitset_vec.size() * temp_bitset_vec2.size();
          contacts_d[i] += temp_bitset_vec_d.size() * temp_bitset_vec_d2.size();
          contacts_d[j] += temp_bitset_vec_d.size() * temp_bitset_vec_d2.size();
          
          // and the running total for the pibd
          pibd[i] += pibd_i;
          pibd[j] += pibd_i;
          pibd_d[i] += pibd_d_i;
          pibd_d[j] += pibd_d_i;
        }
      }
    }
  }
  
  // and now divide accordingly
  for(unsigned int c = 0; c < pop_size; c++) {
    pibd[c] /= (contacts[c] * u_ptr->parameters.g_num_loci);
    pibd_d[c] /= (contacts_d[c] * u_ptr->parameters.g_num_loci);
  }
  
  Rcpp::RawMatrix barcodes = vector_of_bitset_vectors_to_raw_matrix(bitsets, u_ptr->parameters.g_barcode_length);
  
  return(Rcpp::List::create(
      Rcpp::Named("pibd_within") = pibd_within,
      Rcpp::Named("pibd_within_d") = pibd_within_d, 
      Rcpp::Named("pibd") = pibd,
      Rcpp::Named("pibd_d") = pibd_d, 
      Rcpp::Named("age")= ages, 
      Rcpp::Named("state")=states,
      Rcpp::Named("barcodes")=barcodes,
      Rcpp::Named("barcode_states")=barcode_states));
  
}


//' Returns the population's parasite genetics summarised by coi for given sample size and state
//'
//' @param param_list param_list containing statePtr, sample_size, and sample_states
//' @return list of population information
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List population_get_genetics_df_n(Rcpp::List param_list) {
  
  // create universe and do parameter conversions
  Rcpp::XPtr<Universe> u_ptr = Rcpp::as<Rcpp::XPtr<Universe> > (param_list["statePtr"]);
  unsigned int sample_size = Rcpp::as<unsigned int>(param_list["sample_size"]);
  std::vector<unsigned int> sample_states = Rcpp::as<std::vector<unsigned int> >(param_list["sample_states"]);
  
  unsigned int pop_size = u_ptr->population.size();
  
  // handle if we are subsampling
  std::vector<int> pop_sample(pop_size);  
  std::iota(pop_sample.begin(), pop_sample.end(), 0);
  std::vector<unsigned int> all_states = {0, 1, 2, 3, 4, 5};
  
  // if we are subsampling states
  if(all_states != sample_states){
    pop_sample.clear();
    for(unsigned int p = 0 ; p < u_ptr->population.size(); p++) {
      if(std::find(sample_states.begin(), sample_states.end(), static_cast<unsigned int>(u_ptr->population[p].get_m_infection_state())) != sample_states.end()) {
        pop_sample.emplace_back(p);
      }
    }
  }
  
  // if we are choosing a sample
  if(sample_size) {
    pop_sample = sample_without_replacement(pop_sample, sample_size);
    rcpp_out(u_ptr->parameters.g_h_quiet_print, "Sampling");
  }
  pop_size = pop_sample.size();
  
  std::vector< int> coi(pop_size, 0);
  std::vector< int> coi_detected_pcr(pop_size, 0);
  std::vector< int> ages(pop_size, 0);
  std::vector< int> states(pop_size, 0);
  std::vector<std::vector< boost::dynamic_bitset<> > > bitsets(pop_size);
  std::vector<std::vector< int> > barcode_states(pop_size);
  
  // temps
  double q_i;
  double qu_i;
  unsigned int i;
  std::vector<Strain> strains_i;
  std::vector<boost::dynamic_bitset<> > bitsets_pcr_i;
  std::vector<boost::dynamic_bitset<> > bitsets_all_i;
  std::vector< int> barcode_states_i;
  
  
  for (unsigned int el = 0; el < pop_size; el++) {
    
    // who are we sampling
    i = pop_sample[el];
    
    // simple assignments first
    ages[el] = u_ptr->population[i].get_m_person_age();
    states[el] = static_cast< int>(u_ptr->population[i].get_m_infection_state());
    
    // clear temps
    barcode_states_i.clear();
    bitsets_pcr_i.clear();
    bitsets_all_i.clear();
    strains_i.clear();
    
    // reserve space
    strains_i = u_ptr->population[i].get_m_active_strains();
    bitsets_pcr_i.reserve(strains_i.size());
    bitsets_all_i.reserve(strains_i.size());
    barcode_states_i.reserve(strains_i.size());
    barcode_states[el].reserve(strains_i.size());
    bitsets[el].reserve(strains_i.size());
    
    // their detectabilities
    q_i = q_fun(u_ptr->parameters, ages[el], u_ptr->population[i].get_m_ID());
    qu_i = std::pow(q_i,u_ptr->parameters.g_alphaU);
    
    // for coi we want one for all barcodes and then one for just those 
    // that are detected by PCR/Sequenom etc. 
    
    for(auto b : strains_i) {
      if(b.get_m_strain_infection_status() == Strain::SUBPATENT) {
        if(rbernoulli1(qu_i)) {
        bitsets_pcr_i.emplace_back(b.get_m_barcode()); 
        }
      }
      
      bitsets_all_i.emplace_back(b.get_m_barcode()); 
      barcode_states_i.emplace_back(static_cast< int> (b.get_m_strain_infection_status()));
    }
    barcode_states[el] = barcode_states_i;
    bitsets[el] = bitsets_all_i;
    
    // and then the unique size
    if(strains_i.size()){
      
      // cois for pcr
      std::sort(bitsets_pcr_i.begin(), bitsets_pcr_i.end());
      coi_detected_pcr[el] = std::max(1,static_cast<int>(std::unique(bitsets_pcr_i.begin(), bitsets_pcr_i.end()) - bitsets_pcr_i.begin()));

      // similar for all
      std::sort(bitsets_all_i.begin(), bitsets_all_i.end());
      coi[el] = std::unique(bitsets_all_i.begin(), bitsets_all_i.end()) - bitsets_all_i.begin();

    } else {
      
      coi_detected_pcr[el] = 0;
      coi[el] = 0;
      
    }
  }
  
  Rcpp::RawMatrix barcodes = vector_of_bitset_vectors_to_raw_matrix(bitsets, u_ptr->parameters.g_barcode_length);
  
  Rcpp::List res_list = Rcpp::List::create(
    Rcpp::Named("coi", coi),
    Rcpp::Named("coi_detected_pcr",coi_detected_pcr), 
    Rcpp::Named("age",ages), 
    Rcpp::Named("state",states),
    Rcpp::Named("barcodes",barcodes),
    Rcpp::Named("barcode_states",barcode_states));
  
  rcpp_out(u_ptr->parameters.g_h_quiet_print, "post list creation");
  
  return(res_list);
  
}
