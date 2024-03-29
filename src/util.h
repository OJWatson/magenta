// -*-c++-*-
#ifndef UTIL_H
#define UTIL_H

#include <Rcpp.h>
#include <bitset>
#include <R.h>
#include <Rinternals.h>
#include <boost/dynamic_bitset.hpp>
#include "parameters.h"

SEXP bitset_to_sexp(boost::dynamic_bitset<> x, unsigned int n);

SEXP bitset_vector_to_sexp(std::vector<boost::dynamic_bitset<> > x, unsigned int n);

Rcpp::RawVector bitset_vector_to_rawvector(std::vector<boost::dynamic_bitset<> > x, unsigned int n);

Rcpp::RawMatrix bitset_vector_to_raw_matrix(std::vector<boost::dynamic_bitset<> > x, unsigned int n);

Rcpp::RawMatrix vector_of_bitset_vectors_to_raw_matrix(std::vector<std::vector<boost::dynamic_bitset<> > > x, unsigned int n);

boost::dynamic_bitset<> sexp_to_bitset(SEXP x, unsigned int n);

std::vector<std::vector<double> > rcpp_matrix_doubles_to_vec_of_vec(Rcpp::NumericMatrix &mat);

void rcpp_out(bool quiet, std::string messsage);

std::string enum_spatial_convert(Parameters::g_spatial_type_enum v);

#endif
