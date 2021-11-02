#include <cpp11.hpp>
#include <vector>
#include <dust/random/random.hpp>

typedef dust::random::xoroshiro128plus_state rng_state_type;

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

// TODO: all these functions need better names, really
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

  // TODO: this needs implementing:
  // std::vector<double> t_sample =
  //   dust::random::multinom<double>(state, ever_infected_tested, prob);
  std::vector<double>
    t_sample(prob.size(), 1 / static_cast<double>(prob.size()));

  // NOTE: 30 here is an upper bound and is "never going to happen"
  // territory. We expect this to max out about 11 (plus 6 for our
  // negative values)
  std::vector<double> vl_tab(30, 0.0);

  int t_sample_curr = 0;
  int t_sample_seen = 0;
  for (int i = 0; i < ever_infected_tested; ++i) {
    const double a =
      dust::random::normal<double>(state, pars.a_bar, pars.a_sigma);
    const double b =
      dust::random::normal<double>(state, pars.b_bar, pars.b_sigma);
    const double tmax =
      dust::random::normal<double>(state, pars.tmax_bar, pars.tmax_sigma);
    const double log_vlmax =
      dust::random::normal<double>(state, pars.log_vlmax_bar,
                                   pars.log_vlmax_sigma);

    if (t_sample_seen > t_sample[t_sample_curr]) {
      t_sample_curr++;
      t_sample_seen = 0;
    } else {
      t_sample_seen++;
    }

    const auto vl = vl_func(a, b, tmax, t_sample_curr, log_vlmax);
    vl_tab[std::max(vl + 6, 0)]++;
  }

  vl_tab[0] += never_infected_tested;

  double ret = 0.0;
  for (size_t i = 0; i < observed_vl.size(); ++i) {
    ret += observed_vl[i] *
      std::log(vl_tab[i] / static_cast<double>(tested_population));
  }

  return ret;
}

[[cpp11::register]]
cpp11::writable::doubles vl_distribution(cpp11::integers days,
                                         const std::vector<double>& infecteds,
                                         cpp11::list r_observed_vl,
                                         int population,
                                         int tested_population,
                                         cpp11::list r_pars) {
  const vl_parameters pars{cpp11::as_cpp<double>(r_pars["a_bar"]),
                           cpp11::as_cpp<double>(r_pars["a_sigma"]),
                           cpp11::as_cpp<double>(r_pars["b_bar"]),
                           cpp11::as_cpp<double>(r_pars["b_sigma"]),
                           cpp11::as_cpp<double>(r_pars["tmax_bar"]),
                           cpp11::as_cpp<double>(r_pars["tmax_sigma"]),
                           cpp11::as_cpp<double>(r_pars["log_vlmax_bar"]),
                           cpp11::as_cpp<double>(r_pars["log_vlmax_sigma"]) };

  // Cumulative infected over time, just done once. Could accept this
  // as an argument of course but then we rely on it being correct
  // (this way it will always be ok)
  auto cum_infecteds = infecteds;
  for (size_t i = 1; i < infecteds.size(); ++i) {
    cum_infecteds[i] += cum_infecteds[i - 1];
  }

  const int n = days.size();

  // TODO: to do this properly, we should pass the state around
  // (either directly as a vector of raws or as a pointer)
  auto rng = dust::random::prng<rng_state_type>(n, 42);

  // TODO: some validation here that everything is the right size
  // before we access anything. Probably also best to take observed_vl
  // as list.

  std::vector<std::vector<double>> observed_vl;
  for (int i = 0; i < n; ++i) {
    observed_vl.push_back(cpp11::as_cpp<std::vector<double>>(r_observed_vl[i]));
  }

  cpp11::writable::doubles ret(n);
  for (int i = 0; i < n; ++i) {
    ret[i] = vl_distribution_day(rng.state(i),
                                 days[i], infecteds, cum_infecteds,
                                 observed_vl[i],
                                 population, tested_population, pars);
  }

  return ret;
}
