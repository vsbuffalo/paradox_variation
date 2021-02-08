# div_popsize_chains_groups.r -- 
library(tidyverse)
library(scales)
library(latex2exp)
library(HDInterval)


source('plot_funs.r')
source('color_scheme.r')
source('../../R/utilities.r')

opar <- par(no.readonly=TRUE)

# output <- FALSE
output <- TRUE
if (output) 
  pdf('div_popsize_post_phyla.pdf', height=2, width=4)

load('../../data/phyla_diversity_popsize_chains.Rdata')

# we use the separate Stan fit for chordates
# chor_fit <- dt_dps_phyla$fit[[1]]
mols_fit <- dt_dps_phyla$fit[[2]]
arth_fit <- dt_dps_phyla$fit[[3]]

a_phi <- rstan:::extract(arth_fit, pars='phi')$phi
a_beta <- rstan:::extract(arth_fit, pars='beta')$beta

a_sigma_r <- rstan:::extract(arth_fit, pars='sigma_r')$sigma_r 
a_sigma_p <- rstan:::extract(arth_fit, pars='sigma_p')$sigma_p

c_phi <- rstan:::extract(chor_fit, pars='phi')$phi
c_beta <- rstan:::extract(chor_fit, pars='beta')$beta
# no residual error
# c_sigma_r <- rstan:::extract(chor_fit, pars='sigma_r')$sigma_r 
c_sigma_p <- rstan:::extract(chor_fit, pars='sigma_p')$sigma_p

m_phi <- rstan:::extract(mols_fit, pars='phi')$phi
m_beta <- rstan:::extract(mols_fit, pars='beta')$beta
m_sigma_r <- rstan:::extract(mols_fit, pars='sigma_r')$sigma_r 
m_sigma_p <- rstan:::extract(mols_fit, pars='sigma_p')$sigma_p



a_lambda <- a_sigma_p / (a_sigma_p + a_sigma_r)
# c_lambda <- c_sigma_p / (c_sigma_p + c_sigma_r)
m_lambda <- m_sigma_p / (m_sigma_p + m_sigma_r)

layout(matrix(1:9, byrow=TRUE, nrow=3, ncol=3),
       height=c(1, 1, 1.1))
par(mar=c(1.2, 1.4, 1.3, 0.1), oma=c(0.5, 1, 1, 1))

dbs_slope <- dt_dps_phyla %>% 
         ungroup() %>%
         select(slope_lower, slope, slope_upper) %>%
         as.matrix()
dbs_int <- dt_dps_phyla %>% 
         ungroup() %>%
         select(intercept_lower, intercept, intercept_upper) %>%
         as.matrix()
 
tl <- 0.4
post_plot(a_phi, main=latex2exp::TeX('intercept ($\\phi$)'),
          ylab='Arthropoda', bs_ci=dbs_int[3, ], 
          title_line=tl, cex.lab=1, xlim=c(-5, 1))
# text(-4.3, 0.5, 'OLS')
post_plot(a_beta, main=latex2exp::TeX('slope ($\\beta$)'),
          bs_ci=dbs_slope[3, ], xlim=c(-0.4, 0.4))
# text(0.19, 6.5, 'OLS')
post_plot(a_lambda, main=latex2exp::TeX('phylo. signal ($\\lambda$)'),
           xlim=c(0, 1))

post_plot(m_phi, main='',
          ylab='Mollusca', bs_ci=dbs_int[2, ],
          title_line=tl, cex.lab=1, xlim=c(-5, 1))
# text(-4.3, 0.63, 'OLS')
post_plot(m_beta, main='', 
          bs_ci=dbs_slope[2, ], xlim=c(-0.4, 0.4))
# text(0.19, 6.3, 'OLS')
post_plot(m_lambda, main='', xlim=c(0, 1))

par(mar=c(1., 1.4, 1.4, 0.1))
post_plot(c_phi, main='',
          col='seagreen',
          ylab='Chordates', bs_ci=dbs_int[1, ],
          title_line=tl, cex.lab=1, xlim=c(-5, 1))
# text(-4.3, 0.63, 'OLS')
post_plot(c_beta, main='',
          col='seagreen',
          bs_ci=dbs_slope[1, ], xlim=c(-0.4, 0.4))
# text(0.19, 6.3, 'OLS')
post_plot(c_sigma_p^2, 
          col='seagreen',
          main='')
title(main=latex2exp::TeX('phylo. variance ($\\sigma_p^2$)'),
      line=0.3, cex.main=1)

if (output) 
  dev.off()

par(opar)


