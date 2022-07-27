## Direectly taken from Will's code @0e1a6d9, this is the lightly
## cleaned version we were working with.
vl_func <- function(a, b, tmax, t, log_vlmax, k, cap) {
  if(k==2){
    tau = (log_vlmax+tmax)/exp(a)
    vl =  ifelse(t < tau, exp(a)*t-tmax, (1+exp(b)/exp(a))*log_vlmax + exp(b)*tmax/exp(a) - exp(b)*t)
  }

  vl <- floor(-4.87 + 1.418 * log(10^vl))
  vl <- pmin(vl, 30)

}

vl_dist_day2 <- function(day, infecteds, cum_infecteds, population,
                         tested_population,
                         a_params, b_params, tmax_params, log_vlmax_params,
                         negatives, observed, k, cap) {
  proportion_ever_infected <- cum_infecteds[day] / population
  ever_infected_tested <- round(tested_population * proportion_ever_infected, 0)
  never_infected_tested <- tested_population - ever_infected_tested

  positives_tested <- 1000000                                           ######

  # block takes ~ 0.13s
  a_test <- rnorm(ever_infected_tested, a_params[1], a_params[2])
  b_test <- rnorm(ever_infected_tested, b_params[1], b_params[2])
  tmax_test <- rnorm(ever_infected_tested, tmax_params[1], tmax_params[2])
  log_vlmax_test <- rnorm(ever_infected_tested, log_vlmax_params[1],
                          log_vlmax_params[2])


  a_test2 <- rnorm(positives_tested, a_params[1], a_params[2])          #####
  b_test2 <- rnorm(positives_tested, b_params[1], b_params[2])          #####
  tmax_test2 <- rnorm(positives_tested, tmax_params[1], tmax_params[2]) #####
  log_vlmax_test2 <- rnorm(positives_tested, log_vlmax_params[1],       #####
                          log_vlmax_params[2])                          #####


  prob <- infecteds[day:1] / cum_infecteds[day]
  t_sample <- rmultinom(1, ever_infected_tested, prob)
  t_sample_test <- rep(seq_along(t_sample) - 1L, t_sample)

  t_sample2 <- rmultinom(1, positives_tested, prob)                     #####
  t_sample_test2 <- rep(seq_along(t_sample2) - 1L, t_sample2)           #####

  vl <- floor(vl_func(a_test, b_test, tmax_test, t_sample_test,
                      log_vlmax_test, k, cap))

  vl2 <- floor(vl_func(a_test2, b_test2, tmax_test2, t_sample_test2,    #####
                       log_vlmax_test2, k, cap))                        #####

  # block takes ~ 0.1s
  vl_indiv <- data.frame(
    individual = seq_len(ever_infected_tested),
    vl = vl)

  vl_indiv2 <- data.frame(                   #####
    individual = seq_len(positives_tested),  #####
    vl = vl2)                                #####

  ## Several ways of doing this, but this one will be fairly fast
  vl_max <- max(vl)
  vl_tab <- tabulate(pmax(vl, -6) + 7, vl_max + 7)
  vl_tab[[1L]] <- vl_tab[[1L]] + never_infected_tested

  vl_max2 <- max(vl2)                                                                                  #####
  vl_tab2 <- round(tabulate(pmax(vl2, -6) + 7, vl_max2 + 7) * ever_infected_tested/positives_tested,0) #####
  vl_tab2[[1L]] <- vl_tab2[[1L]] + never_infected_tested                                               #####

  len <- max(length(vl_tab), length(observed))
  pad <- function(x, n) {
    if (length(x) < n) {
      c(x, numeric(n - length(x)))
    } else {
      x
    }
  }

  sum(pad(observed, len) * log(pad(vl_tab / tested_population, len)))
  sum(pad(observed, len) * log(pad(vl_tab2 / tested_population, len))) #####

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
                   log_vlmax_bar = 8, log_vlmax_sigma = 2)

  dat$calculate <- function(day, infecteds, observed,
                            population, tested_population, pars) {
      a_params <- c(pars$a_bar, pars$a_sigma)
      b_params <- c(pars$b_bar, pars$b_sigma)
      tmax_params <- c(pars$tmax_bar, pars$tmax_sigma)
      log_vlmax_params <- c(pars$log_vlmax_bar, pars$log_vlmax_sigma)
      cum_infecteds <- cumsum(infecteds)
      vl_dist_day2(day=day,
                   infecteds=infecteds,
                   cum_infecteds=cum_infecteds,
                   population = population,
                   tested_population = population,
                   a_params = a_params,
                   b_params = b_params,
                   tmax_params = tmax_params,
                   log_vlmax_params = log_vlmax_params,
                   negatives = TRUE,
                   observed=observed[[day]]$count,
                   k=2,
                   cap=1)
  }

  dat
}


test_rng_pointer <- function(observed) {
  dust::dust_rng_pointer$new(42, observed$size_full, 1, "xoroshiro128plus")
}
