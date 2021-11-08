## Swapping for ziggurat drops time by about 50% total
##
## Some data dependncy issues are not clear (length of observed vs
## infecteds, etc), is the start day always the first non-empty?

dat <- reference()
dat$cum_infected <- cumsum(dat$infected)
observed <- lapply(dat$observed, function(x)
  as.numeric(x$count) %||% numeric(0))

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

## So dat$cum_infected is a good linear predictor of time:
x <- dat$cum_infected[10:100]
s <- schedule(x, 10)
n <- schedule_n(length(x), 10)
t <- vapply(split(x[s], rep(seq_along(n), n)), sum, numeric(1))
min(t / max(t)) # 99.6%


observed <- tidy_observed(dat$observed)

## 3.411s
system.time(
  sum(vapply(10:100, function(i)
    dat$calculate(i, dat$infected, dat$observed, dat$population,
                  dat$tested_population, dat$pars),
    numeric(1))))

dat$calculate(75, dat$infected, dat$observed, dat$population,
              dat$tested_population, dat$pars)
rng <- rng_init(observed$size, 42)

value <- sapply(10:100, function(i)
  r_calculate_one(i, dat$infected, dat$cum_infected, observed,
                  dat$population, dat$tested_population, dat$pars, rng))

day <- 10:100
value <- sapply(day, function(i)
  system.time(
    r_calculate_one(i, dat$infected, dat$cum_infected, observed,
                    dat$population, dat$tested_population, dat$pars, rng))[[1]])

plot(value ~ day)
plot(value ~ dat$cum_infected[day])

fit <- lm(value ~ dat$cum_infected[day])
fit2 <- lm(tail(value, 51) ~ tail(dat$cum_infected, 51))




tail(n_tasks_threads, n_tasks_ceil * n_threads - n_tasks) <- n_tasks_ceil - 1

k <- rep(m, nt)

sum(k) - length(n)



r_calculate_one(100, dat$infected, dat$cum_infected, observed,
                dat$population, dat$tested_population, dat$pars, rng)
system.time(
r_calculate(10, dat$infected, dat$cum_infected, observed,
            dat$population, dat$tested_population, dat$pars, rng, 1))

## 40 -> 0.264
## 20 -> 0.317
## 10 -> 0.538
system.time(
r_calculate(10, dat$infected, dat$cum_infected, observed,
            dat$population, dat$tested_population, dat$pars, rng, 40))


day <- 10:100
t <- lapply(day, function(i)
  system.time(
    r_calculate_one(i, dat$infected, dat$cum_infected, observed,
                    dat$population, dat$tested_population, dat$pars, rng)))
tt <- sapply(t, "[[", "elapsed")
plot(day, tt)

r_calculate(75L, dat$infected, dat$cum_infected, observed,
            dat$population, dat$tested_population, dat$pars, rng)



bench::mark(
  w = dat$calculate(75, dat$infected, dat$observed, dat$population,
              dat$tested_population, dat$pars),
  r = r_calculate_one(75, dat$infected, dat$cum_infected,
                      dat$observed[[75]]$count,
                      dat$population, dat$tested_population, dat$pars, rng),
  check = FALSE)




vl_distribution(75L, dat$infected, observed, dat$population,
                dat$tested_population, dat$pars)



rng <- rng_init(1, 42)
vl_calculate(75L, dat$infected, dat$cum_infected, dat$observed[[75]]$count,
                   dat$population, dat$tested_population, dat$pars, rng)



list(75L, dat$infected, observed, dat$population,
                dat$tested_population, dat$pars)

## y <- replicate(100,
##                dat$calculate(75, dat$infected, dat$observed, dat$population,
##                              dat$tested_population, dat$pars))

## Here we crash straight away!
