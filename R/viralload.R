##' Compute likelihood
##'
##' @title Compute likelihood
##'
##' @param pars Named list of parameters. Must have elements `a_bar`,
##'   `a_sigma`, `b_bar`, `b_sigma`, `tmax_bar`, `tmax_sigma`,
##'   `log_vlmax_bar`, `log_vlmax_sigma`
##'
##' @param infected Vector of number of simulated infections. Must
##'   have length `observed$size` (this is checked)
##'
##' @param observed An `observed` object, created by `prepare_observed`
##'
##' @param population Integer, size of population
##'
##' @param tested_population Integer, size of tested population
##'
##' @param rng A random number object, created by `rng_init` (this
##'   will swap out soon for a dust function)
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
log_likelihood <- function(pars, infected, observed,
                           population, tested_population, rng,
                           n_threads = 1L, chunk_size = NULL) {
  if (is.null(chunk_size)) {
    chunk_size <- ceiling((observed$size - observed$first) / n_threads)
  }
  if (!inherits(observed, "observed")) {
    stop("Expected an object of class 'observed' for 'observed'")
  }
  if (length(infected) != observed$size) {
    stop(sprintf("Expected observed to have length '%d'", observed$size))
  }
  r_likelihood(pars, infected, observed, population, tested_population, rng,
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
##'   [viralload::log_liklihood]. **DO NOT ALTER THIS OBJECT**
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
  if (length(unique(cutoff)) > 1) {
    stop("Inconsistent negative cut-off")
  }

  cutoff <- -cutoff + 1L
  ret <- list(size = length(len),
              first = which(len > 0)[[1]],
              cutoff = cutoff,
              length = len,
              offset = cumsum(c(0L, len[-length(count)])),
              value = as.integer(unlist(count, FALSE, FALSE)))
  class(ret) <- "observed"
  ret
}

##' Create a streaming RNG object
##'
##' @title Create RNG object
##'
##' @param n_streams Number of independent streams required
##'
##' @param seed The seed (leave `NULL` to use R's seed)
##'
##' @return An opaque pointer to pass through to
##'   [viralload::log_likelihood]
##'
##' @export
rng_init <- function(n_streams, seed = NULL) {
  r_rng_init(n_streams, seed)
}
