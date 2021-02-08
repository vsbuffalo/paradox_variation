# div_popsize_chains.r -- 

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
  pdf('div_popsize_post.pdf', height=1, width=4)

# Stan posteriors, OLS fits, and raw data
load('../../data/diversity_popsize_chains.Rdata')


phi <- rstan:::extract(dps_fit, pars='phi')$phi
beta <- rstan:::extract(dps_fit, pars='beta')$beta
sigma_r <- rstan:::extract(dps_fit, pars='sigma_r')$sigma_r 
sigma_p <- rstan:::extract(dps_fit, pars='sigma_p')$sigma_p

lambda_est(dps_fit)

lambda <- sigma_p / (sigma_p + sigma_r)

layout(matrix(1:3, byrow=TRUE, nrow=1, ncol=3))
par(mar=c(2, 1, 2, 1), oma=rep(0, 4))

post_plot(phi, main=latex2exp::TeX('intercept ($\\phi$)'), 
          bs_ci=dps_bs_cis[1, ])
text(-3.7, 1.23, 'OLS', cex=0.9, xpd=TRUE)
post_plot(beta, main=latex2exp::TeX('slope ($\\beta$)'), 
          bs_ci=dps_bs_cis[2, ])
text(0.125, 20.3, 'OLS', cex=0.9)
post_plot(lambda, main=latex2exp::TeX('phylo. signal ($\\lambda$)'),
           xlim=c(0, 1))

if (output) 
  dev.off()


