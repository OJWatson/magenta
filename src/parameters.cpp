#include "parameters.h"

Parameters::Parameters() {
  
  g_current_time = 1;
  g_calendar_day = 1;
  g_theta = std::vector<double>(365,1.0);
  g_years = 1;
  
  // spatial 
  g_spatial_type = Parameters::NON;
  // metapop spatial stuff
  g_spatial_imported_cotransmission_frequencies = std::vector<int>();
  g_spatial_imported_oocyst_frequencies = std::vector<int>();
  g_spatial_exported_cotransmission_frequencies = std::vector<int>();
  g_spatial_exported_oocyst_frequencies = std::vector<int>();
  g_spatial_total_exported_barcodes = 0;
  g_spatial_total_exported_oocysts = 0;
  g_spatial_total_imported_human_infections = 0;
  g_spatial_total_imported_mosquito_infections = 0;
  g_spatial_imported_human_infection_counter = 0;
  g_spatial_imported_mosquito_infection_counter = 0;
  g_spatial_exported_barcode_counter = 0;
  g_spatial_exported_oocyst_counter = 0;
  g_spatial_exported_barcodes = std::vector<boost::dynamic_bitset<> >();
  g_spatial_exported_oocysts = std::vector<boost::dynamic_bitset<> >();
  g_spatial_imported_barcodes = std::vector<boost::dynamic_bitset<> >();
  g_spatial_imported_oocysts = std::vector<boost::dynamic_bitset<> >();
  
  // more island stuff
  g_percentage_imported_human_infections = 0.0;
  g_percentage_imported_mosquito_infections = 0.0;
  g_total_human_infections = 0;
  g_total_mosquito_infections = 0;
  
  // contransmission and oocyst vectors
  g_cotransmission_frequencies = std::vector<int>();
  g_cotransmission_frequencies_counter = 0;
  g_cotransmission_frequencies_size = 0;
  g_oocyst_frequencies = std::vector<int>();
  g_oocyst_frequencies_counter = 0;
  g_oocyst_frequencies_size = 0;
  
  // demographic parameters
  g_N = 10000;
  g_max_age = 100 * 365;
  g_average_age = 21 * 365;
  // epidemiological parameters
  g_EIR = 200 / 365;
  g_a0 = 8 * 365;
  g_rho = 0.85;
  g_zeta_meanlog = -1.67 / 2.0;
  g_zeta_sdlog = sqrt(1.67);
  // diagnostic parameters
  g_ft = 0.4;
  // entomological parameters
  g_mu0 = 0.132;
  g_mean_mosquito_age = 1/0.132;
  g_beta_gradient = 0.05879106;
  g_beta_intercept = 0.07108894;
  g_ak = 0.30677;
  g_Q0 = 0.92;
  // entomolgical time keeping
  g_mosquito_deficit = 0; 
  g_scourge_today = 0;
  g_mean_mv = 0;
  g_max_mosquito_biting_counter = 1000;
  g_mosquito_next_biting_day_vector = std::vector<int>(g_max_mosquito_biting_counter,3);
  g_mosquito_biting_counter = 0;
  // delays and durations
  g_delay_mos = 10;
  g_mosquito_eip_vector = std::vector<int>(g_max_mosquito_biting_counter, static_cast<int>(g_delay_mos));
  g_mosquito_eip_counter = 0;
  g_mosquito_current_eip = static_cast<int>(g_delay_mos);
  g_delay_gam = 12.5;
  g_dur_E = 12;
  g_dur_T = 5;
  g_dur_D = 5;
  g_dur_U = 110.29;
  g_dur_P = 25;
  g_dur_A = 195;
  g_dur_AU = g_dur_A + g_dur_U;
  // Detection immunity parameters
  g_d1 = 0.160527;
  g_dID = 10 * 365;
  g_ID0 = 1.577533;
  g_kD = 0.476614;
  g_uD = 9.44512;
  g_aD = 8001.99;
  g_fD0 = 0.007055;
  g_gD = 4.8183;
  g_alphaU = 0.185624;
  // Blood stage immunity parameters
  g_b0 = 0.590076;
  g_b1 = 0.5;
  g_dB = 10 * 365;
  g_IB0 = 43.8787;
  g_kB = 2.15506;
  g_uB = 7.19919;
  // Acquired and maternal immunity parameters
  g_phi0 = 0.791666;
  g_phi1 = 0.000737;
  g_dCA = 10950;
  g_IC0 = 18.02366;
  g_kC = 2.36949;
  g_uCA = 6.06349;
  g_PM = 0.774368;
  g_dCM = 67.6952;
  g_mean_maternal_immunity = 0;
  g_sum_maternal_immunity = 0;
  g_total_mums = 0;
  // Severe disease immunity parameters (Griffin et al. 2015 Proc R Soc B 282:20142657, Table S3), using the
  // malariasimulation values. g_fV0, g_aV and g_gammaV are not tabulated there and are the malariasimulation fits
  g_theta0 = 0.0749886;
  g_theta1 = 0.0001191;
  g_dVA = 10950;
  g_IV0 = 1.09629;
  g_kV = 2.00048;
  g_uVA = 11.4321;
  g_fV0 = 0.141195;
  g_aV = 2493.41;
  g_gammaV = 2.91282;
  g_PVM = 0.195768;
  g_dVM = 76.8365;
  g_mean_maternal_severe_immunity = 0;
  g_sum_maternal_severe_immunity = 0;
  // contributions to infectious reservoir by state and age
  g_gamma1 = 1.82425;
  g_cD = 0.0676909;
  g_cT = 0.322 * g_cD;
  g_cU = 0.0062;
  
  // housekeeping parameters
  g_h_quiet_print = false;
  g_h_quiet_test_print = false;
  
  // barcode drug related parameters
  g_resistance_flag = false;
  g_absolute_fitness_cost_flag = false;
  g_number_of_resistance_loci = 0;
  g_resistance_loci = std::vector<unsigned int> {};
  g_artemisinin_loci = {1};
  g_cost_of_resistance = {1.0};
  
  // mutation parameters
  g_mutation_flag = false;
  g_mutation_rate = {0.0000001, 0.0000001};
  g_mutation_treated_modifier = 1.0;
  g_mutations_today = std::vector<unsigned int> {0,0};
  g_mutation_pos_allocator = 0;
  
  // drug related parameters
  g_mft_flag = false;
  g_drug_choice = 0;
  g_number_of_drugs = 0;
  g_partner_drug_ratios = {0.0};
  g_dur_SPC = 10;
  
  // mosquito strain interaction params
  g_vector_adaptation_flag = false;
  g_vector_adaptation_loci = std::vector<unsigned int>{23};
  g_local_oocyst_advantage = 0.5;
  g_gametocyte_sterilisation_flag = false;
  g_gametocyte_sterilisation = 0.5;
  g_oocyst_reduction_by_artemisinin = 0.2;
  
  // non malaria fever parameters
  g_nmf_flag = false;
  g_mean_nmf_frequency = {247.63,232.63,235.94,259.79,298.94,360.32,389.13,446.76,521.43,525.94,475.26,425.41,397.17,361.03}; // probability that a non malarial fever happens
  g_nmf_age_brackets = {-0.1, 365.0, 730.0, 1095.0, 1460.0, 1825.0, 2555.0, 3285.0, 4015.0, 4745.0, 5475.0, 7300.0, 9125.0, 10950.0, 36850.0}; 
  g_prob_of_testing_nmf = 0.5;
    
}  

// barcode global static params
bool Parameters::g_island_imports_plaf_linked_flag = false;
unsigned int Parameters::g_identity_id = 0;
unsigned int Parameters::g_num_loci = 24;
unsigned int Parameters::g_ibd_length = 1;
unsigned int Parameters::g_barcode_length = static_cast<unsigned int>(g_ibd_length * g_num_loci);
std::vector<double> Parameters::g_plaf = std::vector<double> (24,0.5);
std::vector<double> Parameters::g_prob_crossover = std::vector<double> (24,0.5);
Parameters::g_barcode_type_enum Parameters::g_barcode_type = ORDINARY;
