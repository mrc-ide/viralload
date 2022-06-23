##' Compute likelihood
##'
##' @title Compute likelihood
##'
##' @param pars Named list of parameters. Must have elements `a_bar`,
##'   `a_sigma`, `b_bar`, `b_sigma`, `tmax_bar`, `tmax_sigma`,
##'   `log_vlmax_bar`, `log_vlmax_sigma`
##'
##' @param infected Vector of number of simulated infections. Must
##'   have length `observed$size_full` (this is checked)
##'   
##' @param n number of days to truncate infectiousness at
##' 
##' @param k index of the viral load equation to use
##' 
##' @param cap whether to apply a cap to the value of the maximum viral load to avoid NAN errors
##'
##' @param observed An `observed` object, created by `prepare_observed`
##'
##' @param population Integer, size of population
##'
##' @param tested_population Integer, size of tested population
##'
##' @param rng A random number pointer, created by
##'   [dust::dust_rng_pointer]
##'
##' @param n_threads Number of threads to run on, if openmp available.
##'   See [dust::dust_openmp_threads] for guidance
##'
##' @param chunk_size The size of the chunks to run the parallel jobs
##'   in.  The optimal value here will vary depending on your problem
##'   and number of threads. '5' seems like a good guess.
##'
##' @export
##' @return A single numeric log-likelihood value
log_likelihood <- function(pars, n, k, cap, infected, observed,
                           population, tested_population, rng,
                           n_threads = 1L, chunk_size = NULL) {
  if (is.null(chunk_size)) {
    chunk_size <- ceiling(length(observed$day) / n_threads)
  }
  if (!inherits(observed, "observed")) {
    stop("Expected an object of class 'observed' for 'observed'")
  }
  if (length(infected) != observed$size_full) {
    stop(sprintf("Expected observed to have length '%d'", observed$size_full))
  }
  if (length(tested_population) != observed$size_full) {
    if (length(tested_population) == 1L) {
      tested_population <- rep(tested_population, observed$size_full)
    } else {
      stop("Expected 'tested_population' to be a scalar or vector of length ",
           observed$size_full)
    }
  }
  ## Ensure integer storage
  tested_population <- as.integer(tested_population)
  r_likelihood(pars, n, k, cap, infected, observed, population, tested_population, rng,
               n_threads, chunk_size)
}

##' Prepare observed data
##'
##' @title Prepare observed data
##'
##' @param observed A list of `data.frame`s, each with columns `vl`
##'   and `count`
##'
##' @return An object of class `observed` to pass through to
##'   [viralload::log_likelihood]. **DO NOT ALTER THIS OBJECT**
##'
##' @export
prepare_observed <- function(observed) {
  count <- lapply(observed, function(x)
    as.numeric(x$count) %||% numeric(0))

  len <- lengths(count)
  cutoff <- as.integer(vapply(observed[len > 0], function(x) x$vl[[2]], ""))

  cutoff <- NULL
  for (i in which(len > 0)) {
    vl <- observed[[i]]$vl
    if (vl[[1]] != "negative") {
      stop(sprintf("Expected 'negative' for first value of vl in element %d",
                   i))
    }
    if (is.null(cutoff)) {
      cutoff <- as.integer(vl[[2]])
    } else if (as.integer(vl[[2]]) != cutoff) {
      stop(sprintf("Inconsistent first value of vl in element %d, expected %d",
                   i, cutoff))
    }
    stopifnot(identical(
      vl[2:length(vl)],
      as.character(seq(cutoff, length.out = length(vl) - 1L))))
  }

  cutoff <- -cutoff + 1L
  day <- which(len > 0)

  ret <- list(day = day,
              size_full = length(len), # all days of infecteds
              size_data = length(day), # all days we have data for
              cutoff = cutoff,
              length = len,
              offset = cumsum(c(0L, len[-length(len)])),
              value = as.integer(unlist(count, FALSE, FALSE)))
  class(ret) <- "observed"
  ret
}
