// Generated by cpp11: do not edit by hand
// clang-format off


#include "cpp11/declarations.hpp"
#include <R_ext/Visibility.h>

// viralload.cpp
double r_likelihood_one(int r_day, cpp11::list r_pars, int n, int k, int cap, cpp11::doubles r_infecteds, cpp11::list r_viralload, int population, int tested_population, cpp11::sexp r_rng);
extern "C" SEXP _viralload_r_likelihood_one(SEXP r_day, SEXP r_pars, SEXP n, SEXP k, SEXP cap, SEXP r_infecteds, SEXP r_viralload, SEXP population, SEXP tested_population, SEXP r_rng) {
  BEGIN_CPP11
    return cpp11::as_sexp(r_likelihood_one(cpp11::as_cpp<cpp11::decay_t<int>>(r_day), cpp11::as_cpp<cpp11::decay_t<cpp11::list>>(r_pars), cpp11::as_cpp<cpp11::decay_t<int>>(n), cpp11::as_cpp<cpp11::decay_t<int>>(k), cpp11::as_cpp<cpp11::decay_t<int>>(cap), cpp11::as_cpp<cpp11::decay_t<cpp11::doubles>>(r_infecteds), cpp11::as_cpp<cpp11::decay_t<cpp11::list>>(r_viralload), cpp11::as_cpp<cpp11::decay_t<int>>(population), cpp11::as_cpp<cpp11::decay_t<int>>(tested_population), cpp11::as_cpp<cpp11::decay_t<cpp11::sexp>>(r_rng)));
  END_CPP11
}
// viralload.cpp
double r_likelihood(cpp11::list r_pars, int n, int k, int cap, cpp11::doubles r_infecteds, cpp11::list r_viralload, int population, cpp11::integers r_tested_population, cpp11::sexp r_rng, int n_threads, int chunk_size);
extern "C" SEXP _viralload_r_likelihood(SEXP r_pars, SEXP n, SEXP k, SEXP cap, SEXP r_infecteds, SEXP r_viralload, SEXP population, SEXP r_tested_population, SEXP r_rng, SEXP n_threads, SEXP chunk_size) {
  BEGIN_CPP11
    return cpp11::as_sexp(r_likelihood(cpp11::as_cpp<cpp11::decay_t<cpp11::list>>(r_pars), cpp11::as_cpp<cpp11::decay_t<int>>(n), cpp11::as_cpp<cpp11::decay_t<int>>(k), cpp11::as_cpp<cpp11::decay_t<int>>(cap), cpp11::as_cpp<cpp11::decay_t<cpp11::doubles>>(r_infecteds), cpp11::as_cpp<cpp11::decay_t<cpp11::list>>(r_viralload), cpp11::as_cpp<cpp11::decay_t<int>>(population), cpp11::as_cpp<cpp11::decay_t<cpp11::integers>>(r_tested_population), cpp11::as_cpp<cpp11::decay_t<cpp11::sexp>>(r_rng), cpp11::as_cpp<cpp11::decay_t<int>>(n_threads), cpp11::as_cpp<cpp11::decay_t<int>>(chunk_size)));
  END_CPP11
}

extern "C" {
static const R_CallMethodDef CallEntries[] = {
    {"_viralload_r_likelihood",     (DL_FUNC) &_viralload_r_likelihood,     11},
    {"_viralload_r_likelihood_one", (DL_FUNC) &_viralload_r_likelihood_one, 10},
    {NULL, NULL, 0}
};
}

extern "C" attribute_visible void R_init_viralload(DllInfo* dll){
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
  R_forceSymbols(dll, TRUE);
}
