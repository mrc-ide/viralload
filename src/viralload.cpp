#include <cpp11.hpp>
#include <vector>
#include <dust/random/random.hpp>

// TODO: move into dust/interface/random.hpp
#include <cstring>
#include <dust/interface/random.hpp>

// This will move into dust, eventually probably needs to play nice
// with the R-facing RNG as it is literally the same thing.
namespace dust {
namespace random {
template <typename rng_state_type>
SEXP rng_init(int n_threads, cpp11::sexp r_seed) {
  // TODO: we need to pop this into dust/random/r.hpp or similar
  auto seed = dust::interface::as_rng_seed<rng_state_type>(r_seed);
  auto *rng = new prng<rng_state_type>(n_threads, seed);
  auto ret = cpp11::external_pointer<prng<rng_state_type>>(rng);
  return ret;
}

template <typename rng_state_type>
prng<rng_state_type>* rng_get(cpp11::sexp ptr, int n_threads) {
  auto *rng =
    cpp11::as_cpp<cpp11::external_pointer<prng<rng_state_type>>>(ptr).get();
  if (n_threads > 0) {
    if (static_cast<int>(rng->size()) < n_threads) {
      cpp11::stop("Requested a rng with %d threads but only have %d",
                  n_threads, rng->size());
    }
  }
  return rng;
}

}
}

typedef dust::random::xoroshiro128plus_state rng_state_type;

[[cpp11::register]]
cpp11::sexp rng_init(int n_threads, cpp11::sexp seed) {
  return dust::random::rng_init<rng_state_type>(n_threads, seed);
}


struct vl_parameters {
  double a_bar;
  double a_sigma;
  double b_bar;
  double b_sigma;
  double tmax_bar;
  double tmax_sigma;
  double log_vlmax_bar;
  double log_vlmax_sigma;
};

int vl_func(double a, double b, double tmax, double t, double log_vlmax) {
  const auto tau = t - tmax;
  const auto value = std::log10(pow(10, log_vlmax) * (a + b) /
                                (b * std::exp(-a * tau) + a * exp(b * tau)));
  return std::floor(value);
}

double vl_distribution_day(rng_state_type& state,
                           int day,
                           const std::vector<double>& infecteds,
                           const std::vector<double>& cum_infecteds,
                           const std::vector<double>& observed_vl,
                           int population,
                           int tested_population,
                           const vl_parameters& pars) {
  using namespace dust::random;
  // day will come in from R and so be base-1, do the conversion only here
  day--;
  const auto proportion_ever_infected = cum_infecteds[day] / population;
  const auto ever_infected_tested =
    std::round(tested_population * proportion_ever_infected);
  const auto never_infected_tested = tested_population - ever_infected_tested;

  std::vector<double> prob(day + 1);
  for (int i = 0; i <= day; ++i) {
    // TODO: check that day - i is correct here, fairly sure it is.
    // infecteds[day:1] / cum_infecteds[day]
    prob[i] = infecteds[day - i] / cum_infecteds[day];
  }

  std::vector<double> t_sample = multinomial(state, ever_infected_tested, prob);

  // Assume that the first vl_offset cases are the negatives
  const int vl_offset = 6; // fixed for now, could be a parameter
  std::vector<double> vl_tab(observed_vl.size(), 0.0);

  int t_sample_curr = 0;
  int t_sample_seen = 0;
  for (int i = 0; i < ever_infected_tested; ++i) {
    const double a = normal<double>(state, pars.a_bar, pars.a_sigma);
    const double b = normal<double>(state, pars.b_bar, pars.b_sigma);
    const double tmax = normal<double>(state, pars.tmax_bar, pars.tmax_sigma);
    const double log_vlmax = normal<double>(state, pars.log_vlmax_bar,
                                            pars.log_vlmax_sigma);

    if (t_sample_seen > t_sample[t_sample_curr]) {
      t_sample_curr++;
      t_sample_seen = 0;
    } else {
      t_sample_seen++;
    }

    const auto vl = vl_func(a, b, tmax, t_sample_curr, log_vlmax);
    const auto vl_tab_pos = std::max(vl + vl_offset, 0);
    if (static_cast<size_t>(vl_tab_pos) < vl_tab.size()) {
      vl_tab[vl_tab_pos]++;
    }
  }

  vl_tab[0] += never_infected_tested;

  double ret = 0.0;
  for (size_t i = 0; i < observed_vl.size(); ++i) {
    ret += observed_vl[i] *
      std::log(vl_tab[i] / static_cast<double>(tested_population));
  }

  return ret;
}

vl_parameters create_pars(cpp11::list r_pars) {
  return vl_parameters{cpp11::as_cpp<double>(r_pars["a_bar"]),
                         cpp11::as_cpp<double>(r_pars["a_sigma"]),
                         cpp11::as_cpp<double>(r_pars["b_bar"]),
                         cpp11::as_cpp<double>(r_pars["b_sigma"]),
                         cpp11::as_cpp<double>(r_pars["tmax_bar"]),
                         cpp11::as_cpp<double>(r_pars["tmax_sigma"]),
                         cpp11::as_cpp<double>(r_pars["log_vlmax_bar"]),
                         cpp11::as_cpp<double>(r_pars["log_vlmax_sigma"])
                         };
}


double vl_calculate2(rng_state_type& state,
                     int day,
                     const double* infecteds,
                     const double* cum_infecteds,
                     const double* observed_vl,
                     int observed_vl_size,
                     int population,
                     int tested_population,
                     const vl_parameters& pars) {
  using namespace dust::random;
  // day will come in from R and so be base-1, do the conversion only here
  day--;

  // TODO: if we reuse these vectors (infecteds, cum_infecteds) we can
  // compute some of these things ahead of time.
  const auto proportion_ever_infected = cum_infecteds[day] / population;
  const auto ever_infected_tested =
    std::round(tested_population * proportion_ever_infected);
  const auto never_infected_tested = tested_population - ever_infected_tested;

  // TODO: if we support non-normalised prob in multinomial, we can
  // skip this
  std::vector<double> prob(day + 1);
  for (int i = 0; i <= day; ++i) {
    prob[i] = infecteds[day - i] / cum_infecteds[day];
  }

  std::vector<double> t_sample = multinomial(state, ever_infected_tested, prob);

  // Assume that the first vl_offset cases are the negatives
  const int vl_offset = 6; // fixed for now, could be a parameter
  std::vector<double> vl_tab(observed_vl_size, 0.0);

  int t_sample_curr = 0;
  int t_sample_seen = 0;
  for (int i = 0; i < ever_infected_tested; ++i) {
    const double a = normal<double>(state, pars.a_bar, pars.a_sigma);
    const double b = normal<double>(state, pars.b_bar, pars.b_sigma);
    const double tmax = normal<double>(state, pars.tmax_bar, pars.tmax_sigma);
    const double log_vlmax = normal<double>(state, pars.log_vlmax_bar,
                                            pars.log_vlmax_sigma);

    if (t_sample_seen > t_sample[t_sample_curr]) {
      t_sample_curr++;
      t_sample_seen = 0;
    } else {
      t_sample_seen++;
    }

    const auto vl = vl_func(a, b, tmax, t_sample_curr, log_vlmax);
    const auto vl_tab_pos = std::max(vl + vl_offset, 0);
    if (static_cast<size_t>(vl_tab_pos) < vl_tab.size()) {
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

[[cpp11::register]]
double vl_calculate(int day, cpp11::doubles r_infecteds,
                    cpp11::doubles r_cum_infecteds,
                    cpp11::doubles observed,
                    int population,
                    int tested_population,
                    cpp11::list r_pars,
                    cpp11::sexp r_rng) {
  const auto pars = create_pars(r_pars);
  auto rng = dust::random::rng_get<rng_state_type>(r_rng, 1);
  return vl_calculate2(rng->state(0),
                       day,
                       REAL(r_infecteds),
                       REAL(r_cum_infecteds),
                       REAL(observed),
                       observed.size(),
                       population,
                       tested_population,
                       pars);
}

[[cpp11::register]]
cpp11::writable::doubles vl_distribution(cpp11::integers days,
                                         const std::vector<double>& infecteds,
                                         cpp11::list r_observed_vl,
                                         int population,
                                         int tested_population,
                                         cpp11::list r_pars) {
  const auto pars = create_pars(r_pars);

  // Cumulative infected over time, just done once. Could accept this
  // as an argument of course but then we rely on it being correct
  // (this way it will always be ok)
  auto cum_infecteds = infecteds;
  for (size_t i = 1; i < infecteds.size(); ++i) {
    cum_infecteds[i] += cum_infecteds[i - 1];
  }

  const int n = days.size();

  // TODO: to do this properly, we should pass the state around
  // (either directly as a vector of raws or as a pointer). However,
  // we might update to push the parallelisation into the underlying
  // function, which will simplify this considerably.
  //
  // TODO: Add some helpers to create pointers on the R side with some
  // light validation.
  auto rng = dust::random::prng<rng_state_type>(n, 42);

  // TODO: some validation here that everything is the right size
  // before we access anything. Probably also best to take observed_vl
  // as list.
  //
  // TODO: we could do this more efficiently but that may not matter
  // much soon.
  std::vector<std::vector<double>> observed_vl;
  for (int i = 0; i < r_observed_vl.size(); ++i) {
    observed_vl.push_back(cpp11::as_cpp<std::vector<double>>(r_observed_vl[i]));
  }

  cpp11::writable::doubles ret(n);
  for (int i = 0; i < n; ++i) {
    const auto j = days[i] - 1;
    ret[i] = vl_distribution_day(rng.state(i),
                                 days[i], infecteds, cum_infecteds,
                                 observed_vl[j],
                                 population, tested_population, pars);
  }

  return ret;
}
