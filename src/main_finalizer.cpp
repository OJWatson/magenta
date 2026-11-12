//
//  magenta
//  main_finalizer.cpp
//
//  Created: OJ Watson on 06/07/2019
//
//  Distributed under the MIT software licence
//
//  Further details (if any) of this set of functions can be found in the corresponding header file.
//
// ---------------------------------------------------------------------------


#include <iostream>
#include "parameters.h"
#include "probability.h"

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
//' Returns whole model to R in series of nested lists
//'
//' @param param_list parameter list generated with \code{param_list_simulation_finalizer_create}
//' @return list of 1 confirming finalizer has finished
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List Simulation_Finalizer_cpp(Rcpp::List param_list)
{
  
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  // START: R -> C++ CONVERSIONS
  // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  // Create universe pointer from param_list statePtr
  Rcpp::XPtr<Universe> u_ptr = Rcpp::as<Rcpp::XPtr<Universe> > (param_list["statePtr"]);
  
  // prove that C++ code is being run
  rcpp_out(u_ptr->parameters.g_h_quiet_print, "Rcpp function is working!\n");
 
 // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
 // START: MEMORY FREEING
 // ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
 std::vector<Person>().swap(u_ptr->population);
 std::vector<Mosquito>().swap(u_ptr->scourge);
 std::vector<double>().swap(u_ptr->psi_vector);
 std::vector<double>().swap(u_ptr->zeta_vector);
 std::vector<double>().swap(u_ptr->pi_vector);
 u_ptr.release();
 
 return(Rcpp::List::create(Rcpp::Named("Finished") = true));
  
}