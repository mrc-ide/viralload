#include <numeric>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <cpp11.hpp>

#include <dust/random/random.hpp>
#include <dust/r/random.hpp>

// This will move into dust, eventually probably needs to play nice
// with the R-facing RNG as it is literally the same thing.
namespace dust {
namespace random {
namespace r {
template <typename rng_state_type>
SEXP rng_init(int n, cpp11::sexp r_seed) {
  // TODO: we need to pop this into dust/random/r.hpp or similar
  auto seed = dust::r::as_rng_seed<rng_state_type>(r_seed);
  auto *rng = new prng<rng_state_type>(n, seed);
  auto ret = cpp11::external_pointer<prng<rng_state_type>>(rng);
  return ret;
}

template <typename rng_state_type>
prng<rng_state_type>* rng_get(cpp11::sexp ptr, int n) {
  auto *rng =
    cpp11::as_cpp<cpp11::external_pointer<prng<rng_state_type>>>(ptr).get();
  if (n > 0) {
    if (static_cast<int>(rng->size()) < n) {
      cpp11::stop("Requested a rng with %d streams but only have %d",
                  n, rng->size());
    }
  }
  return rng;
}

}
}
}

namespace viralload {

typedef dust::random::xoroshiro128plus_state rng_state_type;

double double_from_list(cpp11::list x, const char * name) {
  return cpp11::as_cpp<double>(x[name]);
}

struct parameters {
  parameters(cpp11::list r_pars) :
    a_bar(double_from_list(r_pars, "a_bar")),
    a_sigma(double_from_list(r_pars, "a_sigma")),
    b_bar(double_from_list(r_pars, "b_bar")),
    b_sigma(double_from_list(r_pars, "b_sigma")),
    tmax_bar(double_from_list(r_pars, "tmax_bar")),
    tmax_sigma(double_from_list(r_pars, "tmax_sigma")),
    log_vlmax_bar(double_from_list(r_pars, "log_vlmax_bar")),
    log_vlmax_sigma(double_from_list(r_pars, "log_vlmax_sigma")) {
  }
  double a_bar;
  double a_sigma;
  double b_bar;
  double b_sigma;
  double tmax_bar;
  double tmax_sigma;
  double log_vlmax_bar;
  double log_vlmax_sigma;
};

struct observed {
  observed(cpp11::list x) :
    size(cpp11::as_cpp<int>(x["size"])),
    length(INTEGER(cpp11::as_cpp<cpp11::integers>(x["length"]))),
    offset(INTEGER(cpp11::as_cpp<cpp11::integers>(x["offset"]))),
    value(INTEGER(cpp11::as_cpp<cpp11::integers>(x["value"]))) {
  }
  const int size;
  const int * length;
  const int * offset;
  const int * value;
};

// TODO: a nicer name here would be good
int vl_func(double a, double b, double tmax, double t, double log_vlmax) {
  const auto tau = t - tmax;
  const auto value = std::log10(pow(10, log_vlmax) * (a + b) /
                                (b * std::exp(-a * tau) + a * exp(b * tau)));
  return std::floor(value);
}

// TODO: if we reuse these vectors (infecteds, cum_infecteds) we can
// compute some of these things ahead of time.
//
// TODO: can pass just 1 cum_infecteds

double calculate_one(const int day,
                     const double* infecteds,
                     const double* cum_infecteds,
                     const observed& viralload,
                     const int population,
                     const int tested_population,
                     const parameters& pars,
                     rng_state_type& state) {
  const auto proportion_ever_infected = cum_infecteds[day] / population;
  const auto ever_infected_tested =
    std::round(tested_population * proportion_ever_infected);
  const auto never_infected_tested = tested_population - ever_infected_tested;

  // NOTE: can be skipped
  std::vector<double> prob(day + 1);
  for (int i = 0; i <= day; ++i) {
    prob[i] = infecteds[day - i] / cum_infecteds[day];
  }

  // Because this will run from openmp, we will crash if this fails
  // for any reason.
  std::vector<int> t_sample(prob.size());
  int * t_sample_data = t_sample.data();
  dust::random::multinomial<double>(state, ever_infected_tested, prob.data(),
                                    prob.size(), t_sample_data);

  const int observed_vl_size = viralload.length[day];
  const int* observed_vl = viralload.value + viralload.offset[day];

  // TODO: lots of integer arithmetic here
  // Assume that the first vl_offset cases are the negatives
  const int vl_offset = 6; // fixed for now, could be a parameter
  std::vector<int> vl_tab(observed_vl_size, 0.0);

  int t_sample_curr = 0;
  int t_sample_seen = 0;
  for (int i = 0; i < ever_infected_tested; ++i) {
    // We'll tidy this up later, see
    // https://github.com/mrc-ide/dust/issues/323
    using dust::random::normal;
    constexpr auto algorithm = dust::random::algorithm::normal::ziggurat;
    const double a =
      normal<double, algorithm>(state, pars.a_bar, pars.a_sigma);
    const double b =
      normal<double, algorithm>(state, pars.b_bar, pars.b_sigma);
    const double tmax =
      normal<double, algorithm>(state, pars.tmax_bar, pars.tmax_sigma);
    const double log_vlmax =
      normal<double, algorithm>(state, pars.log_vlmax_bar,
                                pars.log_vlmax_sigma);

    if (t_sample_seen > t_sample[t_sample_curr]) {
      t_sample_curr++;
      t_sample_seen = 0;
    } else {
      t_sample_seen++;
    }

    const auto vl = vl_func(a, b, tmax, t_sample_curr, log_vlmax);
    const size_t vl_tab_pos = std::max(vl + vl_offset, 0);
    if (vl_tab_pos < vl_tab.size()) {
      vl_tab[vl_tab_pos]++;
    }
  }

  vl_tab[0] += never_infected_tested;

  double ret = 0.0;
  for (int i = 0; i < observed_vl_size; ++i) {
    ret += observed_vl[i] *
      std::log(vl_tab[i] / static_cast<double>(tested_population));
  }

  return ret;
}

double calculate(const int start_day,
                 const double * infecteds,
                 const double * cum_infecteds,
                 const observed& viralload,
                 const int population,
                 const int tested_population,
                 const parameters& pars,
                 dust::random::prng<rng_state_type>* rng,
                 int n_threads,
                 int chunk_size) {
  const int len = viralload.size - start_day + 1;
  double ret = 0.0;
#ifdef _OPENMP
#pragma omp parallel for schedule(static, chunk_size) num_threads(n_threads) reduction(+:ret)
#endif
  for (int i = 0; i < len; ++i) {
    auto& state = rng->state(i);
    const int day = viralload.size - i - 1;
    ret += calculate_one(day, infecteds, cum_infecteds,
                         viralload, population, tested_population,
                         pars, state);
  }
  return ret;
}

}

[[cpp11::register]]
cpp11::sexp rng_init(int n_threads, cpp11::sexp seed) {
  return dust::random::r::rng_init<viralload::rng_state_type>(n_threads, seed);
}

// Interface to calculate just one day; we'll use this for debugging mostly
[[cpp11::register]]
double r_calculate_one(int r_day,
                       cpp11::doubles r_infecteds,
                       cpp11::doubles r_cum_infecteds,
                       cpp11::list r_viralload,
                       int population,
                       int tested_population,
                       cpp11::list r_pars,
                       cpp11::sexp r_rng) {
  const int day = r_day - 1;
  const double * infecteds = REAL(r_infecteds);
  const double * cum_infecteds = REAL(r_cum_infecteds);
  const viralload::observed viralload(r_viralload);
  const viralload::parameters pars(r_pars);
  auto rng = dust::random::r::rng_get<viralload::rng_state_type>(r_rng, 1);
  return viralload::calculate_one(day, infecteds, cum_infecteds,
                                  viralload, population, tested_population,
                                  pars, rng->state(0));
}

[[cpp11::register]]
double r_calculate(int r_start_day,
                   cpp11::doubles r_infecteds,
                   cpp11::doubles r_cum_infecteds,
                   cpp11::list r_viralload,
                   int population,
                   int tested_population,
                   cpp11::list r_pars,
                   cpp11::sexp r_rng,
                   int n_threads,
                   int chunk_size) {
  // assert observed.size is same as (cum_)infecteds size
  // assert start_day within range
  const int start_day = r_start_day;
  const double * infecteds = REAL(r_infecteds);
  const double * cum_infecteds = REAL(r_cum_infecteds);
  const viralload::observed viralload(r_viralload);
  const viralload::parameters pars(r_pars);
  auto rng = dust::random::r::rng_get<viralload::rng_state_type>(r_rng, 1);

  return viralload::calculate(start_day, infecteds, cum_infecteds,
                              viralload, population, tested_population,
                              pars, rng, n_threads, chunk_size);
}
