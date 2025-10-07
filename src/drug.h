//
//  magenta
//  drug.h
//
//  Created: OJ on 13/01/2017
//	Most recent edits: OJ on 10/03/2017
//
//  Distributed under the MIT software licence - see Notes.c file for details
//
//  Class describing a parasite strain.
//
// ---------------------------------------------------------------------------


#ifndef DRUG_H
#define DRUG_H

#include <vector>
#include <queue>
#include <bitset>
#include <random>
#include <boost/dynamic_bitset.hpp>
#include "probability.h"
#include "Rcpp.h"

class Drug {
  
public:
  
  
private:
  
  std::vector<double> m_lpf; // prob of lpf for each bitset combination related to this drug moda
  std::vector<unsigned int> m_barcode_positions; // which positions in the barcode correspond to this drug
  std::vector<unsigned int> m_prophylactic_positions; // which positions in the barcode correspond to the prophylactic
  double m_dur_P;					// duration of prophylaxis
  double m_dur_SPC;					// duration of slow parasite clearance
  double m_hill_n; // Hill function parameter n for curve detailing pretective efficacy
  double m_hill_kA; // Hill function parameter kA for curve detailing pretective efficacy
  double m_hill_res_n; // Hill function parameter n for curve detailing pretective efficacy when challenged by resistant parasite
  double m_hill_res_kA; // Hill function parameter kA for curve detailing pretective efficacy when challenged by resistant parasite
  
public:
  
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // CONSTRUCTORS
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  // Class constructor using required parameters
  Drug(std::vector<double> lpf, 
       std::vector<unsigned int> barcode_positions, 
       std::vector<unsigned int> prophylactic_positions, 
       double dur_P, 
       double dur_SPC,
       double hill_n,
       double hill_kA,
       double hill_res_n,
       double hill_res_kA) : 
  
  m_lpf(lpf), 
  m_barcode_positions(barcode_positions),
  m_prophylactic_positions(prophylactic_positions),
  m_dur_P(dur_P),
  m_dur_SPC(dur_SPC),
  m_hill_n(hill_n),
  m_hill_kA(hill_kA),
  m_hill_res_n(hill_res_n),
  m_hill_res_kA(hill_res_kA)
  {};

  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // GETTERS
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  // Get lpf
  std::vector<double> get_m_lpf() const { return(m_lpf); }
  
  // Get barcode positions
  std::vector<unsigned int> get_m_barcode_positions() const { return(m_barcode_positions); }
  
  // Get prophylactic positions
  std::vector<unsigned int> get_m_prophylactic_positions() const { return(m_prophylactic_positions); }
  
  // Get drug duration of prophylaxis
  double get_m_dur_P() const { return(m_dur_P); }		
  
  // Getduration of slow parasite clearance
  double get_m_dur_SPC() const { return(m_dur_SPC); }		
  
  // Get Hill function parameter n for curve detailing pretective efficacy
  double get_m_hill_n() const { return(m_hill_n); }	
  
  // Get Hill function parameter kA for curve detailing pretective efficacy
  double get_m_hill_kA() const { return(m_hill_kA); }	
  
  // Get Hill function parameter n for curve detailing pretective efficacy when challenged by resistant parasite
  double get_m_hill_res_n() const { return(m_hill_res_n); }	
  
  // Get Hill function parameter kA for curve detailing pretective efficacy when challenged by resistant parasite
  double get_m_hill_res_kA() const { return(m_hill_res_kA); }	
  
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // Setters
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  // Set lpf
  void set_m_lpf(std::vector<double> x) { m_lpf = x; }
  
  // Set barcode positions
  void set_m_barcode_positions(std::vector<unsigned int> x) { m_barcode_positions = x; }
  
  // Set prophylactic positions
  void set_m_prophylactic_positions(std::vector<unsigned int> x) { m_prophylactic_positions = x; }
  
  // Set drug duration of prophylaxis
  void set_m_dur_P(double x) { m_dur_P = x; }
  
  // Set  duration of slow parasite clearance
  void set_m_dur_SPC(double x) { m_dur_SPC = x; }
  
  // Set Hill function parameter n for curve detailing pretective efficacy
  void set_m_hill_n(double x) { m_hill_n = x; }	
  
  // Set Hill function parameter kA for curve detailing pretective efficacy
  void set_m_hill_kA(double x) { m_hill_kA = x; }	
  
  // Set Hill function parameter n for curve detailing pretective efficacy when challenged by resistant parasite
  void set_m_hill_res_n(double x) { m_hill_res_n = x; }	
  
  // Set Hill function parameter kA for curve detailing pretective efficacy when challenged by resistant parasite
  void set_m_hill_res_kA(double x) { m_hill_res_kA = x; }	
  
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // Extra
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  // Get prob of lpf
  double get_prob_of_lpf_x(unsigned int x) const { return(m_lpf[x]); }
  
  // Get prob of lpf at position
  double get_prob_of_lpf_barcode(boost::dynamic_bitset<> &x) const { 
    
    unsigned long mask = 1;
    unsigned long result = 0;
    for (unsigned int start_bit = 0 ; start_bit < m_barcode_positions.size(); start_bit++) {
      if (x[m_barcode_positions[start_bit]]) {
        result |= mask;
      }
      mask <<= 1;
    }
    return (m_lpf[result]);
    
    }
  
  // Does the provided barcode have resistant loci
  bool resistance_to_drug(boost::dynamic_bitset<> &x) const { 
    
    // start as false and then and assign
    bool resistant = false;
    
    // what is the strain's lpf and def(ault)
    double def = get_prob_of_lpf_x(0);
    double lpf = get_prob_of_lpf_barcode(x);
    
    // if it is lower then check if this is because it has prophylactic resistance positions
    if (lpf < def) {
    
    for (unsigned int p : m_prophylactic_positions) {
      resistant = resistant || x[p];
    }
    
    }
    
    return(resistant);
    
  }
  
  // Get early reinfection
  bool early_reinfection(boost::dynamic_bitset<> &x,
                         unsigned int current_time,
                         unsigned int day_of_change,
                         unsigned int day_treated) const { 
    
    // is the infecting strain resistant at prophylactic positions
    bool resistant = resistance_to_drug(x);
    
    // what is their drug concentration
    double drug_conc = R::dexp(current_time - day_treated, 
                               1.0/(day_of_change - day_treated), 
                               false);
    
    // normalise so that drug conc is 1 at t = 0
    drug_conc = drug_conc / R::dexp(0.0, 1.0/(day_of_change - day_treated), false);
    
    // if it was resistant then check for early reinfection using resistant hill parameters
    if (resistant) {
      
      return(rbernoulli1(hill_function(drug_conc, m_hill_res_n, m_hill_res_kA)));
      
    } else {
      
      return(rbernoulli1(hill_function(drug_conc, m_hill_n, m_hill_kA)));
      
    }
    
  }
  
  
  
  //
  
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // RcppList Conversion
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  // Class converter
  Rcpp::List drug_to_rcpp_list() {
    return(
      Rcpp::List::create(
        Rcpp::Named("m_lpf")=m_lpf, 
        Rcpp::Named("m_barcode_positions")=m_barcode_positions,
        Rcpp::Named("m_prophylactic_positions")=m_prophylactic_positions,
        Rcpp::Named("m_dur_P")=m_dur_P,
        Rcpp::Named("m_dur_SPC")=m_dur_SPC,
        Rcpp::Named("m_hill_n")=m_hill_n,
        Rcpp::Named("m_hill_kA")=m_hill_kA,
        Rcpp::Named("m_hill_res_n")=m_hill_res_n,
        Rcpp::Named("m_hill_res_kA")=m_hill_res_kA
      )
    );
  }
  
  // Class constructor from RcppList
  Drug(Rcpp::List list) : 
    
    m_lpf(Rcpp::as<std::vector<double> >(list["m_lpf"])), 
    m_barcode_positions(Rcpp::as<std::vector<unsigned int> >(list["m_barcode_positions"])), 
    m_prophylactic_positions(Rcpp::as<std::vector<unsigned int> >(list["m_prophylactic_positions"])), 
    m_dur_P(Rcpp::as<double>(list["m_dur_P"])),
    m_dur_SPC(Rcpp::as<double>(list["m_dur_SPC"])),
    m_hill_n(Rcpp::as<double>(list["m_hill_n"])),
    m_hill_kA(Rcpp::as<double>(list["m_hill_kA"])),
    m_hill_res_n(Rcpp::as<double>(list["m_hill_res_n"])),
    m_hill_res_kA(Rcpp::as<double>(list["m_hill_res_kA"]))
    
  {};
  
  
  
};

#endif
