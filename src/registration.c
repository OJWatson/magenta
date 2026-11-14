#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void odin_itn_irs_initmod_desolve(void *);
extern void odin_itn_irs_output_dde(void *);
extern void odin_itn_irs_rhs_dde(void *);
extern void odin_itn_irs_rhs_desolve(void *);

/* .Call calls */
extern SEXP _magenta_population_get_genetics_df_n(SEXP);
extern SEXP _magenta_population_get_genetics_ibd_df_n(SEXP);
extern SEXP _magenta_Simulation_Finalizer_cpp(SEXP);
extern SEXP _magenta_Simulation_Get_cpp(SEXP);
extern SEXP _magenta_Simulation_Init_cpp(SEXP);
extern SEXP _magenta_Simulation_Saved_Init_cpp(SEXP);
extern SEXP _magenta_Simulation_Update_cpp(SEXP);
extern SEXP _magenta_test_barcode_from_PLAF(SEXP, SEXP);
extern SEXP _magenta_test_bitset_serialisation(SEXP, SEXP);
extern SEXP _magenta_test_dependent_barcode_from_PLAF(SEXP, SEXP);
extern SEXP _magenta_test_generate_next_ibd(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _magenta_test_ibd_conversion(SEXP, SEXP, SEXP, SEXP);
extern SEXP _magenta_test_recombinant_with_ibd(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _magenta_test_relative_contribution(SEXP, SEXP, SEXP, SEXP);
extern SEXP odin_itn_irs_contents(SEXP);
extern SEXP odin_itn_irs_create(SEXP);
extern SEXP odin_itn_irs_initial_conditions(SEXP, SEXP);
extern SEXP odin_itn_irs_metadata(SEXP);
extern SEXP odin_itn_irs_rhs_r(SEXP, SEXP, SEXP);
extern SEXP odin_itn_irs_set_initial(SEXP, SEXP, SEXP, SEXP);
extern SEXP odin_itn_irs_set_user(SEXP, SEXP);

static const R_CMethodDef CEntries[] = {
    {"odin_itn_irs_initmod_desolve", (DL_FUNC) &odin_itn_irs_initmod_desolve, 1},
    {"odin_itn_irs_output_dde",      (DL_FUNC) &odin_itn_irs_output_dde,      1},
    {"odin_itn_irs_rhs_dde",         (DL_FUNC) &odin_itn_irs_rhs_dde,         1},
    {"odin_itn_irs_rhs_desolve",     (DL_FUNC) &odin_itn_irs_rhs_desolve,     1},
    {NULL, NULL, 0}
};

static const R_CallMethodDef CallEntries[] = {
    {"_magenta_population_get_genetics_df_n",     (DL_FUNC) &_magenta_population_get_genetics_df_n,     1},
    {"_magenta_population_get_genetics_ibd_df_n", (DL_FUNC) &_magenta_population_get_genetics_ibd_df_n, 1},
    {"_magenta_Simulation_Finalizer_cpp",         (DL_FUNC) &_magenta_Simulation_Finalizer_cpp,         1},
    {"_magenta_Simulation_Get_cpp",               (DL_FUNC) &_magenta_Simulation_Get_cpp,               1},
    {"_magenta_Simulation_Init_cpp",              (DL_FUNC) &_magenta_Simulation_Init_cpp,              1},
    {"_magenta_Simulation_Saved_Init_cpp",        (DL_FUNC) &_magenta_Simulation_Saved_Init_cpp,        1},
    {"_magenta_Simulation_Update_cpp",            (DL_FUNC) &_magenta_Simulation_Update_cpp,            1},
    {"_magenta_test_barcode_from_PLAF",           (DL_FUNC) &_magenta_test_barcode_from_PLAF,           2},
    {"_magenta_test_bitset_serialisation",        (DL_FUNC) &_magenta_test_bitset_serialisation,        2},
    {"_magenta_test_dependent_barcode_from_PLAF", (DL_FUNC) &_magenta_test_dependent_barcode_from_PLAF, 2},
    {"_magenta_test_generate_next_ibd",           (DL_FUNC) &_magenta_test_generate_next_ibd,           5},
    {"_magenta_test_ibd_conversion",              (DL_FUNC) &_magenta_test_ibd_conversion,              4},
    {"_magenta_test_recombinant_with_ibd",        (DL_FUNC) &_magenta_test_recombinant_with_ibd,        6},
    {"_magenta_test_relative_contribution",       (DL_FUNC) &_magenta_test_relative_contribution,       4},
    {"odin_itn_irs_contents",                     (DL_FUNC) &odin_itn_irs_contents,                     1},
    {"odin_itn_irs_create",                       (DL_FUNC) &odin_itn_irs_create,                       1},
    {"odin_itn_irs_initial_conditions",           (DL_FUNC) &odin_itn_irs_initial_conditions,           2},
    {"odin_itn_irs_metadata",                     (DL_FUNC) &odin_itn_irs_metadata,                     1},
    {"odin_itn_irs_rhs_r",                        (DL_FUNC) &odin_itn_irs_rhs_r,                        3},
    {"odin_itn_irs_set_initial",                  (DL_FUNC) &odin_itn_irs_set_initial,                  4},
    {"odin_itn_irs_set_user",                     (DL_FUNC) &odin_itn_irs_set_user,                     2},
    {NULL, NULL, 0}
};

void R_init_magenta(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
