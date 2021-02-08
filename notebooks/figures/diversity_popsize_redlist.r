# diversity_popsize_redlist.r -- 

library(tidyverse)
library(plotrix)
library(ggrepel)
library(wesanderson)
library(mixtools)
library(scales)
library(sp)
library(grid)
library(latex2exp)

source('plot_funs.r')

opar <- par(no.readonly=TRUE)


# for stan fit
load('../../data/diversity_popsize_chains.Rdata')
load('../../data/main_datasets.Rdata')
phi <- mean(rstan:::extract(dps_fit, pars='phi')$phi)
beta <- mean(rstan:::extract(dps_fit, pars='beta')$beta)

# alias this -- less typing
d <- da_dps

iucn_levels <- c("NA/DD", "LC", "NT", "VU", "EN", "CR")
cols <- c("gray70", "#006666", "aquamarine3", "#cc9900", "#cc6633", "#cc3333")
d$redlist_cat <- factor(d$redlist_cat, levels=iucn_levels)

# what to include
# x <- d$pred_log10_N  # the simplistic unregularized version
x <- d$log10_popsize
y <- d$log10_diversity

xlims <- c(4, 18)
plot(x, y)

output <- TRUE
# output <- FALSE

if (output)
  pdf("diversity_popsize_redlist.pdf", width=7, height=7/1.5)

par(mar=c(4, 4, 3, 4), oma=c(0.5, 0.5, 0.5, 0.5))

d$redlist_cat[is.na(d$redlist_cat)] <- "NA/DD"
d$redlist_cat[d$redlist_cat == 'DD'] <- "NA/DD" # there's only two, just convert
#iucn_cols <- setNames(c('gray80', wes_palette("Cavalcanti1")), iucn_levels) 
iucn_cols <- setNames(cols, iucn_levels) 

plot(x, y, col=iucn_cols[d$redlist_cat], type='n', axes=FALSE,
     ylab='', xlab='', ylim=c(-4.1, 0), xlim=xlims) 

adjust <- 0
points(x, y, pch=21, cex=1.3, lwd=0.4,
       bg=iucn_cols[d$redlist_cat], col='white')
xseq <- seq(4, 18, 2)
axis(1, xseq, 
     labels=latex2exp::TeX(sprintf("$10^{%d}$", xseq)), 
     padj = -1,
     tck=-0.01,
     cex.axis=0.8, line=0.5)
yseq <- seq(-4, 0)
axis(2, yseq, 
     labels=latex2exp::TeX(sprintf("$10^{%d}$", yseq)), 
     las=1,
     tck=-0.01, hadj=0.5,
     cex.axis=0.8, line=0.5)

y <- -0
m <- 0.12
x <- 18
n <- 0.25
for (i in seq_along(iucn_levels)) {
  status <- iucn_levels[[i]]
  # pop size range
  u <- d$log10_popsize[d$redlist_cat == status]
  rng <- quantile(u, probs=c(0.25, 0.75))
  segments(rng[1], y + i*m, rng[2], y+i*m, col=iucn_cols[status], lwd=3, lend=2, xpd=TRUE)
  points(mean(u), y + i*m, pch=19, col=iucn_cols[status], xpd=TRUE)
  rng <- range(u)
  segments(rng[1], y + i*m, rng[2], y+i*m, col=iucn_cols[status], lwd=1, lend=2, xpd=TRUE)

  # diversity range
  u <- d$log10_diversity[d$redlist_cat == status]
  rng <- quantile(u, probs=c(0.25, 0.75))
  segments(x + i*n, rng[1], x+i*n, rng[2], col=iucn_cols[status], lwd=3, lend=2, xpd=TRUE)
  points(x + i*n, mean(u), pch=19, col=iucn_cols[status], xpd=TRUE)
  rng <- range(u)
  segments(x + i*n, rng[1], x+i*n, rng[2], col=iucn_cols[status], lwd=1, lend=2, xpd=TRUE)
}

# new dataset with factors
drl <- da_dps
drl$redlist_cat[is.na(drl$redlist_cat)] <- "NA/DD"
# ordered in terms of endangeredness
iucn_levels <- c("NA/DD", "LC", "NT", "VU", "EN", "CR")
drl$redlist_cat <- factor(drl$redlist_cat, levels=iucn_levels)





# rseq <- yseq - log10(4e-8)
# axis(4, yseq, 
#      labels=latex2exp::TeX(sprintf("$10^{%.2f}$", rseq)), 
#      las=1,
#      tck=-0.01, hadj=0.5,
#      cex.axis=0.8, line=0.5)
title(ylab="   pairwise diversity (differences per bp)", 
      line=2.8, cex.lab=1.2)
title(xlab=latex2exp::TeX("approximate population size"), line=2.5, cex.lab=1.2)


iucn_cols_s <- iucn_cols[iucn_levels]
legend(8, -4, names(iucn_cols_s), fill=iucn_cols_s,
       x.intersp=1,
       bty='n', border=0, cex=0.6, ncol=6)

if (output) dev.off()



par(opar)
