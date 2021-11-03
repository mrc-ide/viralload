## pkgload::load_all("~/Documents/Projects/ncov/dust", export_all = FALSE)
pkgload::load_all()

dat <- reference()
observed <- lapply(dat$observed, function(x)
  as.numeric(x$count) %||% numeric(0))
vl_distribution(75L, dat$infected, observed, dat$population,
                  dat$tested_population, dat$pars)

dat$cum_infected <- cumsum(dat$infected)

rng <- rng_init(1, 42)
vl_calculate(75L, dat$infected, dat$cum_infected, dat$observed[[75]]$count,
                   dat$population, dat$tested_population, dat$pars, rng)


bench::mark(
  a = vl_distribution(75L, dat$infected, observed, dat$population,
                  dat$tested_population, dat$pars),
  b = dat$calculate(75, dat$infected, dat$observed, dat$population,
                    dat$tested_population, dat$pars),
  c = vl_calculate(75L, dat$infected, dat$cum_infected, dat$observed[[75]]$count,
                   dat$population, dat$tested_population, dat$pars, rng),
  check = FALSE)

list(75L, dat$infected, observed, dat$population,
                dat$tested_population, dat$pars)

## y <- replicate(100,
##                dat$calculate(75, dat$infected, dat$observed, dat$population,
##                              dat$tested_population, dat$pars))

## Here we crash straight away!
