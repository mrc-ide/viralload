## NOTE: it looks to me that I've probably got a small amount of bias
## in the normal calculation here when using the ziggurat algorithm;
## if that concerns you then I'd suggest switching to use box_muller
## (search for dust::random::algorithm in the C++ code). The effect is
## small but can be seen when plotting a distribution. However, your
## code from previously is occassionally producing NaN values so I'd
## suspect that there could be a bug lurking there too.
test_that("calculation agrees with the R version for one day", {
  dat <- reference()
  dat$tested_population <- 1e6
  observed <- prepare_observed(dat$observed)
  rng <- test_rng_pointer(observed)
  day <- 30L
  y <- replicate(
    30,
    r_likelihood_one(day, dat$pars, n=30, k=2, cap=30, dat$infected, observed,
                     dat$population, dat$tested_population, rng))
  set.seed(1)
  cmp <- replicate(
    30,
    dat$calculate(day=day, infecteds=dat$infected, observed=dat$observed,
                  population=dat$population, tested_population=dat$tested_population,
                  pars=dat$pars))

  expect_gt(t.test(y, cmp)$p.value, 0.05)
})


## Final calculation is deterministically the sum of the parts
test_that("calculation is sum over days", {
  dat <- reference()
  observed <- prepare_observed(dat$observed)

  rng1 <- test_rng_pointer(observed)
  rng2 <- test_rng_pointer(observed)

  ll1 <- log_likelihood(dat$pars, n=30, k=2, cap=30, dat$infected, observed,
                        dat$population, dat$tested_population,
                        rng1)
  ll2 <- vapply(observed$day, r_likelihood_one, numeric(1),
                dat$pars, n=30, k=2, cap=30, dat$infected, observed,
                dat$population, dat$tested_population, rng2)
  expect_equal(ll1, sum(ll2))
})


test_that("calculation is sum over days when testing varies", {
  dat <- reference()
  observed <- prepare_observed(dat$observed)

  rng1 <- test_rng_pointer(observed)
  rng2 <- test_rng_pointer(observed)

  tested_population <- rpois(observed$size_full, 1e6)

  ll1 <- log_likelihood(dat$pars, n=30, k=2, cap=30, dat$infected, observed,
                        dat$population, tested_population,
                        rng1)
  ll2 <- vapply(10:100, function(day)
    r_likelihood_one(day, dat$pars, n=30, k=2, cap=30, dat$infected, observed,
                     dat$population, tested_population[day], rng2),
    numeric(1))

  expect_equal(ll1, sum(ll2))
})


test_that("Check that infected is the correct size", {
  dat <- reference()
  observed <- prepare_observed(dat$observed)
  rng <- test_rng_pointer(observed)
  expect_error(
    log_likelihood(dat$pars, n=30, k=2, cap=30, dat$infected[-1], observed,
                   dat$population, dat$tested_population,
                   rng),
    "Expected observed to have length '100'")
})


test_that("Check that we are given 'observed' object", {
  dat <- reference()
  observed <- prepare_observed(dat$observed)
  rng <- test_rng_pointer(observed)
  expect_error(
    log_likelihood(dat$pars, dat$infected, NULL,
                   dat$population, dat$tested_population,
                   rng),
    "Expected an object of class 'observed' for 'observed'")
})


test_that("Check that cuttoff is consistent", {
  dat <- reference()
  obs <- dat$observed
  obs[[52]] <- obs[[52]][-2, ]
  expect_error(
    prepare_observed(obs),
    "Inconsistent first value of vl in element 52, expected -5")
  obs[[13]] <- obs[[13]][-1, ]
  expect_error(
    prepare_observed(obs),
    "Expected 'negative' for first value of vl in element 13")
})


test_that("tested_population must be of suitable length", {
  dat <- reference()
  observed <- prepare_observed(dat$observed)

  rng <- test_rng_pointer(observed)
  tested_population <- rep(1e6, observed$size_full - 1)

  expect_error(
    log_likelihood(dat$pars, n=30, k=2, cap=1, dat$infected, observed,
                   dat$population, tested_population,
                   rng),
    "Expected 'tested_population' to be a scalar or vector of length 100")
})


test_that("calculation is sum over days", {
  dat <- reference()

  i <- sample.int(100, 10)
  observed_empty <- dat$observed
  observed_empty[i] <- list(NULL)
  observed_empty <- prepare_observed(observed_empty)

  observed_full <- prepare_observed(dat$observed)

  ## TODO: what is the expectation here?
  rng1 <- test_rng_pointer(observed_empty)
  ll1 <- log_likelihood(dat$pars, n=30, k=2, cap=1, dat$infected, observed_empty,
                       dat$population, dat$tested_population,
                       rng1)

  rng2 <- test_rng_pointer(observed_full)
  ll2 <- vapply(observed_empty$day, r_likelihood_one, n=30, k=2, cap=1, numeric(1),
                dat$pars, dat$infected, observed_full,
                dat$population, dat$tested_population, rng2)

  expect_equal(ll1, sum(ll2))
})
