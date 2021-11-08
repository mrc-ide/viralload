tidy_observed <- function(d) {
  d <- lapply(dat$observed, function(x)
    as.numeric(x$count) %||% numeric(0))
  len <- lengths(d)
  list(size = length(len),
       length = len,
       offset = cumsum(c(0L, len[-length(d)])),
       value = unlist(d, FALSE, FALSE)) # int?
}

schedule <- function(x, n_threads) {
  if (n_threads == 1L) {
    return(seq_along(x))
  }

  n_tasks <- schedule_n(length(x), n_threads)

  n <- integer(n_threads)
  time <- numeric(n_threads)
  tasks <- lapply(n_tasks, integer)

  idx <- order(x, decreasing = TRUE)
  for (i in seq_along(idx)) {
    task <- idx[[i]]
    pos <- n < n_tasks
    j <- which(pos)[which.min(time[pos])]
    k <- n[[j]] <- n[[j]] + 1L
    tasks[[j]][[k]] <- task
    time[[j]] <- time[[j]] + x[[task]]
  }

  unlist(tasks)
}


schedule_n <- function(n_tasks, n_threads) {
  min <- floor(n_tasks / n_threads)
  ret <- rep(min, n_threads)
  ret[seq_len(n_tasks - min * n_threads)] <- min + 1L
  ret
}
