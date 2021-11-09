// Generated by cpp11: do not edit by hand
// clang-format off


#include "cpp11/declarations.hpp"
#include <R_ext/Visibility.h>

// viralload.cpp
cpp11::sexp rng_init(int n_threads, cpp11::sexp seed);
extern "C" SEXP _viralload_rng_init(SEXP n_threads, SEXP seed) {
  BEGIN_CPP11
    return cpp11::as_sexp(rng_init(cpp11::as_cpp<cpp11::decay_t<int>>(n_threads), cpp11::as_cpp<cpp11::decay_t<cpp11::sexp>>(seed)));
  END_CPP11
}
// viralload.cpp
double r_calculate_one(int r_day, cpp11::doubles r_infecteds, cpp11::doubles r_cum_infecteds, cpp11::list r_viralload, int population, int tested_population, cpp11::list r_pars, cpp11::sexp r_rng);
extern "C" SEXP _viralload_r_calculate_one(SEXP r_day, SEXP r_infecteds, SEXP r_cum_infecteds, SEXP r_viralload, SEXP population, SEXP tested_population, SEXP r_pars, SEXP r_rng) {
  BEGIN_CPP11
    return cpp11::as_sexp(r_calculate_one(cpp11::as_cpp<cpp11::decay_t<int>>(r_day), cpp11::as_cpp<cpp11::decay_t<cpp11::doubles>>(r_infecteds), cpp11::as_cpp<cpp11::decay_t<cpp11::doubles>>(r_cum_infecteds), cpp11::as_cpp<cpp11::decay_t<cpp11::list>>(r_viralload), cpp11::as_cpp<cpp11::decay_t<int>>(population), cpp11::as_cpp<cpp11::decay_t<int>>(tested_population), cpp11::as_cpp<cpp11::decay_t<cpp11::list>>(r_pars), cpp11::as_cpp<cpp11::decay_t<cpp11::sexp>>(r_rng)));
  END_CPP11
}
// viralload.cpp
double r_calculate(int r_start_day, cpp11::doubles r_infecteds, cpp11::doubles r_cum_infecteds, cpp11::list r_viralload, int population, int tested_population, cpp11::list r_pars, cpp11::sexp r_rng, int n_threads, int chunk_size);
extern "C" SEXP _viralload_r_calculate(SEXP r_start_day, SEXP r_infecteds, SEXP r_cum_infecteds, SEXP r_viralload, SEXP population, SEXP tested_population, SEXP r_pars, SEXP r_rng, SEXP n_threads, SEXP chunk_size) {
  BEGIN_CPP11
    return cpp11::as_sexp(r_calculate(cpp11::as_cpp<cpp11::decay_t<int>>(r_start_day), cpp11::as_cpp<cpp11::decay_t<cpp11::doubles>>(r_infecteds), cpp11::as_cpp<cpp11::decay_t<cpp11::doubles>>(r_cum_infecteds), cpp11::as_cpp<cpp11::decay_t<cpp11::list>>(r_viralload), cpp11::as_cpp<cpp11::decay_t<int>>(population), cpp11::as_cpp<cpp11::decay_t<int>>(tested_population), cpp11::as_cpp<cpp11::decay_t<cpp11::list>>(r_pars), cpp11::as_cpp<cpp11::decay_t<cpp11::sexp>>(r_rng), cpp11::as_cpp<cpp11::decay_t<int>>(n_threads), cpp11::as_cpp<cpp11::decay_t<int>>(chunk_size)));
  END_CPP11
}

extern "C" {
static const R_CallMethodDef CallEntries[] = {
    {"_viralload_r_calculate",     (DL_FUNC) &_viralload_r_calculate,     10},
    {"_viralload_r_calculate_one", (DL_FUNC) &_viralload_r_calculate_one,  8},
    {"_viralload_rng_init",        (DL_FUNC) &_viralload_rng_init,         2},
    {NULL, NULL, 0}
};
}

extern "C" attribute_visible void R_init_viralload(DllInfo* dll){
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
  R_forceSymbols(dll, TRUE);
}
