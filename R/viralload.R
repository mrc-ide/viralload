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


prepare_observed <- function(d) {
  count <- lapply(d, function(x)
    as.numeric(x$count) %||% numeric(0))

  len <- lengths(count)
  cutoff <- as.integer(vapply(d[len > 0], function(x) x$vl[[2]], ""))

  cutoff <- NULL
  for (i in which(len > 0)) {
    vl <- d[[i]]$vl
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
