// Generated by cpp11: do not edit by hand
// clang-format off


#include "cpp11/declarations.hpp"
#include <R_ext/Visibility.h>

// viralload.cpp
double vl_calculate(int day, cpp11::doubles r_infecteds, cpp11::doubles r_cum_infecteds, cpp11::doubles observed, int population, int tested_population, cpp11::list r_pars);
extern "C" SEXP _viralload_vl_calculate(SEXP day, SEXP r_infecteds, SEXP r_cum_infecteds, SEXP observed, SEXP population, SEXP tested_population, SEXP r_pars) {
  BEGIN_CPP11
    return cpp11::as_sexp(vl_calculate(cpp11::as_cpp<cpp11::decay_t<int>>(day), cpp11::as_cpp<cpp11::decay_t<cpp11::doubles>>(r_infecteds), cpp11::as_cpp<cpp11::decay_t<cpp11::doubles>>(r_cum_infecteds), cpp11::as_cpp<cpp11::decay_t<cpp11::doubles>>(observed), cpp11::as_cpp<cpp11::decay_t<int>>(population), cpp11::as_cpp<cpp11::decay_t<int>>(tested_population), cpp11::as_cpp<cpp11::decay_t<cpp11::list>>(r_pars)));
  END_CPP11
}
// viralload.cpp
cpp11::writable::doubles vl_distribution(cpp11::integers days, const std::vector<double>& infecteds, cpp11::list r_observed_vl, int population, int tested_population, cpp11::list r_pars);
extern "C" SEXP _viralload_vl_distribution(SEXP days, SEXP infecteds, SEXP r_observed_vl, SEXP population, SEXP tested_population, SEXP r_pars) {
  BEGIN_CPP11
    return cpp11::as_sexp(vl_distribution(cpp11::as_cpp<cpp11::decay_t<cpp11::integers>>(days), cpp11::as_cpp<cpp11::decay_t<const std::vector<double>&>>(infecteds), cpp11::as_cpp<cpp11::decay_t<cpp11::list>>(r_observed_vl), cpp11::as_cpp<cpp11::decay_t<int>>(population), cpp11::as_cpp<cpp11::decay_t<int>>(tested_population), cpp11::as_cpp<cpp11::decay_t<cpp11::list>>(r_pars)));
  END_CPP11
}

extern "C" {
static const R_CallMethodDef CallEntries[] = {
    {"_viralload_vl_calculate",    (DL_FUNC) &_viralload_vl_calculate,    7},
    {"_viralload_vl_distribution", (DL_FUNC) &_viralload_vl_distribution, 6},
    {NULL, NULL, 0}
};
}

extern "C" attribute_visible void R_init_viralload(DllInfo* dll){
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
  R_forceSymbols(dll, TRUE);
}
