## NOTE: it looks to me that I've probably got a small amount of bias
## in the normal calculation here when using the ziggurat algorithm;
## if that concerns you then I'd suggest switching to use box_muller
## (search for dust::random::algorithm in the C++ code). The effect is
## small but can be seen when plotting a distribution. However, your
## code from previously is occassionally producing NaN values so I'd
## suspect that there could be a bug lurking there too.
test_that("calculation agrees with the R version for one day", {
  dat <- reference()
  observed <- prepare_observed(dat$observed)
  rng <- test_rng_pointer(observed)
  day <- 30L
  y <- replicate(
    30,
    r_likelihood_one(day, dat$pars, dat$infected, observed,
                     dat$population, dat$tested_population, rng))
  set.seed(1)
  cmp <- replicate(
    30,
    dat$calculate(day, dat$infected, dat$observed,
                  dat$population, dat$tested_population,
                  dat$pars))

  expect_gt(t.test(y, cmp)$p.value, 0.05)
})


## Final calculation is deterministically the sum of the parts
test_that("calculation is sum over days", {
  dat <- reference()
  observed <- prepare_observed(dat$observed)

  rng1 <- test_rng_pointer(observed)
  rng2 <- test_rng_pointer(observed)

  ll1 <- log_likelihood(dat$pars, dat$infected, observed,
                        dat$population, dat$tested_population,
                        rng1)
  ll2 <- vapply(10:100, r_likelihood_one, numeric(1),
                dat$pars, dat$infected, observed,
                dat$population, dat$tested_population, rng2)
  expect_equal(ll1, sum(ll2))
})


test_that("calculation is sum over days when testing varies", {
  dat <- reference()
  observed <- prepare_observed(dat$observed)

  rng1 <- test_rng_pointer(observed)
  rng2 <- test_rng_pointer(observed)

  tested_population <- rpois(observed$size, 1e6)

  ll1 <- log_likelihood(dat$pars, dat$infected, observed,
                        dat$population, tested_population,
                        rng1)
  ll2 <- vapply(10:100, function(day)
    r_likelihood_one(day, dat$pars, dat$infected, observed,
                     dat$population, tested_population[day], rng2),
    numeric(1))

  expect_equal(ll1, sum(ll2))
})


test_that("Check that infected is the correct size", {
  dat <- reference()
  observed <- prepare_observed(dat$observed)
  rng <- test_rng_pointer(observed)
  expect_error(
    log_likelihood(dat$pars, dat$infected[-1], observed,
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
  tested_population <- rep(1e6, observed$size - 1)

  expect_error(
    log_likelihood(dat$pars, dat$infected, observed,
                   dat$population, tested_population,
                   rng),
    "Expected 'tested_population' to be a scalar or vector of length 100")
})
