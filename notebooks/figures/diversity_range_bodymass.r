library(tidyverse)
source('color_scheme.r')
source('../../R/utilities.r')

load('../../data/main_datasets.Rdata')
load('../../data/diversity_range_bodymass_chains.Rdata')

d <- da_dps

output <- TRUE
# output <- FALSE

if (output)
  pdf("diversity_range_bodymass.pdf", width=7, height=3)


par(mar=c(0.5, 0.5, 1, 0.5), oma=rep(0.4, 4))
opar <- par(no.readonly=TRUE)
layout(matrix(c(1, 2, 3, 4, 5, 6, 7,
                8, 8, 8, 9, 10, 10, 10), 
              byrow=TRUE, ncol=7, nrow=2),
      height=c(0.3, 1), width=c(0.3, 0.3, 0.3, 0.1, 0.3, 0.3, 0.3))


phi <- rstan:::extract(dbm_fit, pars='phi')$phi
beta <- rstan:::extract(dbm_fit, pars='beta')$beta
sigma_r <- rstan:::extract(dbm_fit, pars='sigma_r')$sigma_r 
sigma_p <- rstan:::extract(dbm_fit, pars='sigma_p')$sigma_p

lambda <- sigma_p / (sigma_p + sigma_r)

post_plot(phi, main=latex2exp::TeX('intercept ($\\phi$)'))
# text(-3.7, 1.23, 'OLS', cex=0.9)
post_plot(beta, main=latex2exp::TeX('slope ($\\beta$)'))
# text(0.140, 20.3, 'OLS', cex=0.9)
post_plot(lambda, main=latex2exp::TeX('phylo. signal ($\\lambda$)'),
           xlim=c(0, 1))

par(mar=rep(0, 4))
plot.new()


par(mar=c(0.5, 0.5, 1, 0.5))
phi <- rstan:::extract(drn_fit, pars='phi')$phi
beta <- rstan:::extract(drn_fit, pars='beta')$beta
sigma_r <- rstan:::extract(drn_fit, pars='sigma_r')$sigma_r 
sigma_p <- rstan:::extract(drn_fit, pars='sigma_p')$sigma_p

lambda <- sigma_p / (sigma_p + sigma_r)

post_plot(phi, main=latex2exp::TeX('intercept ($\\phi$)'))
# text(-3.7, 1.23, 'OLS', cex=0.9)
post_plot(beta, main=latex2exp::TeX('slope ($\\beta$)'))
# text(0.140, 20.3, 'OLS', cex=0.9)
post_plot(lambda, main=latex2exp::TeX('phylo. signal ($\\lambda$)'),
           xlim=c(0, 1))


par(mar=c(3, 4, 1, 0))

y <- d$log10_diversity
x <- d$log10_body_mass
plot(x, y, col=phyla_cols[d$phylum], type='n', 
     axes=FALSE,
     ylab='', xlab='', ylim=c(-4.1, 0)) 

points(x, y, pch=21, cex=1.3, lwd=0.4,
       bg=phyla_cols[d$phylum], col='white')
xseq <- seq(-4, 6, 2)
yseq <- (-4):0
axis(1, xseq, 
     line=0.3,
     padj = -0.9,
     cex.axis=0.8,
     tck=-0.02,
     labels=latex2exp::TeX(sprintf("$10^{%d}$", xseq)))

axis(2, yseq, 
     las=1, 
     cex.axis=0.8,
     tck=-0.02, hadj=0.65,
     labels=latex2exp::TeX(sprintf("$10^{%d}$", yseq)),
     line=0.5)
mu <- colMeans(rstan:::extract(dbm_fit, 'mu_rep')$mu_rep)
z <- log10_body_mass_seq
lines(z, mu, lwd=1.2, lty=2, col='cornflowerblue')
lines(z, coef(dbm_ols)[1] + coef(dbm_ols)[2]*z, lwd=1.2, lty=2, col='gray42')

title(ylab="pairwise diversity", line=2.8, cex.lab=1.2)
title(xlab=latex2exp::TeX("body mass (g)"), line=2, cex.lab=1.2)

par(mar=rep(0, 4))
plot.new()


par(mar=c(3, 4, 1, 0))
y <- d$log10_diversity
x <- d$log10_range
plot(x, y, col=phyla_cols[d$phylum], type='n', 
     axes=FALSE,
     ylab='', xlab='', ylim=c(-4.1, 0)) 
xseq <- seq(3, 8)
axis(1, xseq, 
     line=0.3,
     padj = -0.9,
     cex.axis=0.8,
     tck=-0.02,
     labels=latex2exp::TeX(sprintf("$10^{%d}$", xseq)))
axis(2, yseq, 
     las=1, 
     cex.axis=0.8,
     tck=-0.02, hadj=0.65,
     labels=latex2exp::TeX(sprintf("$10^{%d}$", yseq)),
     line=0.5)

points(x, y, pch=21, cex=1.3, lwd=0.4,
       bg=phyla_cols[d$phylum], col='white')

mu <- colMeans(rstan:::extract(drn_fit, 'mu_rep')$mu_rep)
z <- log10_range_seq
lines(z, mu, lwd=1.2, lty=2, col='cornflowerblue')
lines(z, coef(drn_ols)[1] + coef(drn_ols)[2]*z, lwd=1.2, lty=2, col='gray42')

title(ylab="pairwise diversity", line=2.8, cex.lab=1.2)
title(xlab=latex2exp::TeX("range (km^2)"), line=2, cex.lab=1.2)

legend(2.6, 0, names(phyla_cols), fill=phyla_cols,
       bty='n', border=0, cex=0.7, ncol=3)



if (output)
  dev.off()

