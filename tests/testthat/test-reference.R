reference <- function() {
  vl_dist_list <- vl_dist_list_generator(infecteds_df = infecteds)

  list(
    calculate = function(day, infecteds, observed_vl, population,
                         tested_population, pars) {
      a_params <- c(pars$a_bar, pars$a_sigma)
      b_params <- c(pars$b_bar, pars$b_sigma)
      tmax_params <- c(pars$tmax_bar, pars$tmax_sigma)
      log_vlmax_params <- c(pars$log_vlmax_bar, pars$log_vlmax_sigma)
      cum_infecteds <- cumsum(infecteds)
      vl_dist_day2(day,
                   sim_infecteds$infecteds,
                   sim_infecteds$cum_infecteds,
                   population = population,
                   tested_population = population,
                   a_params = a_params,
                   b_params = b_params,
                   tmax_params = tmax_params,
                   log_vlmax_params = log_vlmax_params,
                   negatives = TRUE,
                   observed_vl)
    })
}
