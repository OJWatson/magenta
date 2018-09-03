// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// Simulation_Get_cpp
Rcpp::List Simulation_Get_cpp(Rcpp::List paramList);
RcppExport SEXP _MAGENTA_Simulation_Get_cpp(SEXP paramListSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type paramList(paramListSEXP);
    rcpp_result_gen = Rcpp::wrap(Simulation_Get_cpp(paramList));
    return rcpp_result_gen;
END_RCPP
}
// Simulation_Init_cpp
Rcpp::List Simulation_Init_cpp(Rcpp::List paramList);
RcppExport SEXP _MAGENTA_Simulation_Init_cpp(SEXP paramListSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type paramList(paramListSEXP);
    rcpp_result_gen = Rcpp::wrap(Simulation_Init_cpp(paramList));
    return rcpp_result_gen;
END_RCPP
}
// Simulation_Saved_Init_cpp
Rcpp::List Simulation_Saved_Init_cpp(Rcpp::List paramList);
RcppExport SEXP _MAGENTA_Simulation_Saved_Init_cpp(SEXP paramListSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type paramList(paramListSEXP);
    rcpp_result_gen = Rcpp::wrap(Simulation_Saved_Init_cpp(paramList));
    return rcpp_result_gen;
END_RCPP
}
// Simulation_Update_cpp
Rcpp::List Simulation_Update_cpp(Rcpp::List paramList);
RcppExport SEXP _MAGENTA_Simulation_Update_cpp(SEXP paramListSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type paramList(paramListSEXP);
    rcpp_result_gen = Rcpp::wrap(Simulation_Update_cpp(paramList));
    return rcpp_result_gen;
END_RCPP
}
// test_barcode_from_PLAF
SEXP test_barcode_from_PLAF(Rcpp::NumericVector plaf, unsigned int n);
RcppExport SEXP _MAGENTA_test_barcode_from_PLAF(SEXP plafSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type plaf(plafSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(test_barcode_from_PLAF(plaf, n));
    return rcpp_result_gen;
END_RCPP
}
// test_recombinant_with_ibd
SEXP test_recombinant_with_ibd(SEXP barcode_1, SEXP barcode_2, unsigned int bl, unsigned int nl, unsigned int ib, Rcpp::NumericVector pc);
RcppExport SEXP _MAGENTA_test_recombinant_with_ibd(SEXP barcode_1SEXP, SEXP barcode_2SEXP, SEXP blSEXP, SEXP nlSEXP, SEXP ibSEXP, SEXP pcSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type barcode_1(barcode_1SEXP);
    Rcpp::traits::input_parameter< SEXP >::type barcode_2(barcode_2SEXP);
    Rcpp::traits::input_parameter< unsigned int >::type bl(blSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type nl(nlSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type ib(ibSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type pc(pcSEXP);
    rcpp_result_gen = Rcpp::wrap(test_recombinant_with_ibd(barcode_1, barcode_2, bl, nl, ib, pc));
    return rcpp_result_gen;
END_RCPP
}
// test_generate_next_ibd
SEXP test_generate_next_ibd(unsigned int bl, unsigned int nl, unsigned int ib, Rcpp::NumericVector pc, unsigned int id);
RcppExport SEXP _MAGENTA_test_generate_next_ibd(SEXP blSEXP, SEXP nlSEXP, SEXP ibSEXP, SEXP pcSEXP, SEXP idSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< unsigned int >::type bl(blSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type nl(nlSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type ib(ibSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type pc(pcSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type id(idSEXP);
    rcpp_result_gen = Rcpp::wrap(test_generate_next_ibd(bl, nl, ib, pc, id));
    return rcpp_result_gen;
END_RCPP
}
// test_ibd_conversion
Rcpp::List test_ibd_conversion(SEXP barcode, unsigned int bl, unsigned int nl, unsigned int ib);
RcppExport SEXP _MAGENTA_test_ibd_conversion(SEXP barcodeSEXP, SEXP blSEXP, SEXP nlSEXP, SEXP ibSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type barcode(barcodeSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type bl(blSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type nl(nlSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type ib(ibSEXP);
    rcpp_result_gen = Rcpp::wrap(test_ibd_conversion(barcode, bl, nl, ib));
    return rcpp_result_gen;
END_RCPP
}
// population_get_genetics_ibd_df_n
Rcpp::List population_get_genetics_ibd_df_n(Rcpp::List paramList);
RcppExport SEXP _MAGENTA_population_get_genetics_ibd_df_n(SEXP paramListSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type paramList(paramListSEXP);
    rcpp_result_gen = Rcpp::wrap(population_get_genetics_ibd_df_n(paramList));
    return rcpp_result_gen;
END_RCPP
}
// population_get_genetics_df_n
Rcpp::List population_get_genetics_df_n(Rcpp::List paramList);
RcppExport SEXP _MAGENTA_population_get_genetics_df_n(SEXP paramListSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type paramList(paramListSEXP);
    rcpp_result_gen = Rcpp::wrap(population_get_genetics_df_n(paramList));
    return rcpp_result_gen;
END_RCPP
}
// test_bitset_serialisation
SEXP test_bitset_serialisation(SEXP x, unsigned int n);
RcppExport SEXP _MAGENTA_test_bitset_serialisation(SEXP xSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type x(xSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(test_bitset_serialisation(x, n));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_MAGENTA_Simulation_Get_cpp", (DL_FUNC) &_MAGENTA_Simulation_Get_cpp, 1},
    {"_MAGENTA_Simulation_Init_cpp", (DL_FUNC) &_MAGENTA_Simulation_Init_cpp, 1},
    {"_MAGENTA_Simulation_Saved_Init_cpp", (DL_FUNC) &_MAGENTA_Simulation_Saved_Init_cpp, 1},
    {"_MAGENTA_Simulation_Update_cpp", (DL_FUNC) &_MAGENTA_Simulation_Update_cpp, 1},
    {"_MAGENTA_test_barcode_from_PLAF", (DL_FUNC) &_MAGENTA_test_barcode_from_PLAF, 2},
    {"_MAGENTA_test_recombinant_with_ibd", (DL_FUNC) &_MAGENTA_test_recombinant_with_ibd, 6},
    {"_MAGENTA_test_generate_next_ibd", (DL_FUNC) &_MAGENTA_test_generate_next_ibd, 5},
    {"_MAGENTA_test_ibd_conversion", (DL_FUNC) &_MAGENTA_test_ibd_conversion, 4},
    {"_MAGENTA_population_get_genetics_ibd_df_n", (DL_FUNC) &_MAGENTA_population_get_genetics_ibd_df_n, 1},
    {"_MAGENTA_population_get_genetics_df_n", (DL_FUNC) &_MAGENTA_population_get_genetics_df_n, 1},
    {"_MAGENTA_test_bitset_serialisation", (DL_FUNC) &_MAGENTA_test_bitset_serialisation, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_MAGENTA(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
