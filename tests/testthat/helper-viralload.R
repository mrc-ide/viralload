## Direectly taken from Will's code @0e1a6d9, this is the lightly
## cleaned version we were working with.
vl_func <- function(a, b, tmax, t, log_vlmax) {
  log10(10^log_vlmax * (a+b) / (b * exp(-a*(t-tmax)) + a*exp(b*(t-tmax))))
}

vl_dist_day2 <- function(day, infecteds, cum_infecteds, population,
                         tested_population,
                         a_params, b_params, tmax_params, log_vlmax_params,
                         negatives, observed) {
  proportion_ever_infected <- cum_infecteds[day] / population
  ever_infected_tested <- round(tested_population * proportion_ever_infected, 0)
  never_infected_tested <- tested_population - ever_infected_tested

  # block takes ~ 0.13s
  a_test <- rnorm(ever_infected_tested, a_params[1], a_params[2])
  b_test <- rnorm(ever_infected_tested, b_params[1], b_params[2])
  tmax_test <- rnorm(ever_infected_tested, tmax_params[1], tmax_params[2])
  log_vlmax_test <- rnorm(ever_infected_tested, log_vlmax_params[1],
                          log_vlmax_params[2])

  prob <- infecteds[day:1] / cum_infecteds[day]
  t_sample <- rmultinom(1, ever_infected_tested, prob)
  t_sample_test <- rep(seq_along(t_sample) - 1L, t_sample)

  vl <- floor(vl_func(a_test, b_test, tmax_test, t_sample_test,
                      log_vlmax_test))
  # block takes ~ 0.1s
  vl_indiv <- data.frame(
    individual = seq_len(ever_infected_tested),
    vl = vl)

  ## Several ways of doing this, but this one will be fairly fast
  vl_max <- max(vl)
  vl_tab <- tabulate(pmax(vl, -6) + 7, vl_max + 7)
  vl_tab[[1L]] <- vl_tab[[1L]] + never_infected_tested

  len <- max(length(vl_tab), length(observed))
  pad <- function(x, n) {
    if (length(x) < n) {
      c(x, numeric(n - length(x)))
    } else {
      x
    }
  }

  sum(pad(observed, len) * log(pad(vl_tab / tested_population, len)))
}


reference <- function() {
  ## To generate the data from will's code @0e1a6d9
  ## saveRDS(list(infected = infecteds,
  ##              observed = vl_dist_list),
  ##         "example.rds", version = 2L)

  dat <- readRDS("example.rds")
  dat$population <- 1e6
  dat$tested_population <- 1e6
  dat$pars <- list(a_bar = 3.0, a_sigma = 0.4,
                   b_bar = 1.0, b_sigma = 0.1,
                   tmax_bar = 7, tmax_sigma = 1,
                   log_vlmax_bar = 7, log_vlmax_sigma = 1)

  dat$calculate <- function(day, infecteds, observed,
                            population, tested_population, pars) {
      a_params <- c(pars$a_bar, pars$a_sigma)
      b_params <- c(pars$b_bar, pars$b_sigma)
      tmax_params <- c(pars$tmax_bar, pars$tmax_sigma)
      log_vlmax_params <- c(pars$log_vlmax_bar, pars$log_vlmax_sigma)
      cum_infecteds <- cumsum(infecteds)
      vl_dist_day2(day,
                   infecteds,
                   cum_infecteds,
                   population = population,
                   tested_population = population,
                   a_params = a_params,
                   b_params = b_params,
                   tmax_params = tmax_params,
                   log_vlmax_params = log_vlmax_params,
                   negatives = TRUE,
                   observed[[day]]$count)
  }

  dat
}
