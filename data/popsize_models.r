library(tidyverse)
library(rstan)

# x prefixes are so this can be included in combined_dataset.r

# My combined dataset, animal subset
xda <- d %>% #read_tsv('../data/combined_data.tsv') %>%
       filter(kingdom == 'Animalia') %>%
       filter_at(vars(size, range), ~ is.finite(.)) %>%
       # put a minimum range of 1km2
       mutate(range = ifelse(range == 0, 1, range))

# Romiguier et al (2014) data
xdr <- dr %>%
        #read_tsv('../data/romiguier_et_al_2014_updated.tsv') %>%
        filter_at(vars(size, bodymass_g), ~ !is.na(.))

# Damuth 1987 data
xdd <- dd #read_tsv('../data/damuth1987.tsv')

plot(log10(bodymass_g) ~ log10(size), data=xdr)
plot(log10(density) ~ log10(mass_g), data=xdd)

dat <- list( # Romiguier et al data
             N_r = nrow(xdr),
             body_length_r=xdr$size,
             body_mass_r=xdr$bodymass_g,
             # Damuth 1987 data
             N_d = nrow(xdd),
             density_d = xdd$density,
             body_mass_d = xdd$mass_g,
             # main dataset
             N = nrow(xda),
             body_length = xda$size,
             range_meas = xda$range)


species <- xda$species

CACHED_POPSIZE_CHAINS_FILE <- "popsize_chains.Rdata"
FORCE <- FALSE
if (FORCE || !file.exists(CACHED_POPSIZE_CHAINS_FILE)) {
  bayes_fit <- stan('../stan/pred_popsize_missing_centered.stan', data=dat, 
                    cores=4,
                    iter=10000)
  save(bayes_fit, file=CACHED_POPSIZE_CHAINS_FILE)
} else {
  load(CACHED_POPSIZE_CHAINS_FILE)
}

pars <- c('alpha', 'beta', 'sigma_mass', 'gamma', 'theta', 'sigma_density')
print(bayes_fit, pars=pars)
# print(fit, pars='pred_density')

popsize <- rstan:::extract(bayes_fit, pars='log10_popsize')$log10_popsize
mass <- rstan:::extract(bayes_fit, pars='log10_body_mass')$log10_body_mass
density <- rstan:::extract(bayes_fit, pars='log10_density')$log10_density

# sampled from normal_rng using parameters
mass_pred <- rstan:::extract(bayes_fit, pars='log10_body_mass_pred')$log10_body_mass_pred
density_pred <- rstan:::extract(bayes_fit, pars='log10_density_pred')$log10_density_pred

ndraws <- ncol(popsize)

species_mean <- function(x) lapply(1:ncol(x), function(i) x[, i])

post_popsizes <- tibble(species = species, 
                       log10_popsize = species_mean(popsize),
                       log10_body_mass = species_mean(mass),
                       log10_density = species_mean(density))

save(post_popsizes, file='post_pred_draws.Rdata')

bayes_popsizes <- tibble(species = species, 
                         log10_body_mass = colMeans(mass),
                         log10_body_mass_pred = colMeans(mass_pred),
                         log10_popsize = colMeans(popsize),
                         log10_density = colMeans(density),
                         log10_density_pred = colMeans(density_pred))

# traceplot(bayes_fit, pars=pars)


dn <- colMeans(extract(bayes_fit, pars='log10_density')$log10_density)
bm <- colMeans(extract(bayes_fit, pars='log10_body_mass')$log10_body_mass)
ps <- colMeans(extract(bayes_fit, pars='log10_popsize')$log10_popsize)
psr <- extract(bayes_fit, pars='log10_popsize_range')$log10_popsize_range

