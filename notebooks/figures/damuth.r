library(tidyverse)
library(scales)

opar <- par(no.readonly=TRUE)

pdf("damuth.pdf", width=5, height=5)
par(mar=c(5, 5, 2, 1))

load('../../data/popsize_chains.Rdata')

pars <- rstan:::extract(bayes_fit, pars=c('theta', 'gamma'))
mean_slope <- mean(pars$theta)
mean_intercept <- mean(pars$gamma)

d <- read_tsv('../../data/damuth1987.tsv')  %>%
       mutate(group=factor(group)) %>%
       mutate(log10_mass_g = log10(mass_g),
              log10_density = log10(density))

ngroups <- length(unique(d$group))
group_cols <- hue_pal()(ngroups)

plot(d$log10_mass_g, d$log10_density, pch=19, 
     type='n', axes=FALSE, xlab='',
     ylab='')
points(d$log10_mass_g, d$log10_density, pch=21, 
       cex=1, lwd=0.4, bg=group_cols[as.factor(d$group)], col='white')
abline(a=mean_intercept, b=mean_slope, lty=2, col='gray52')
xseq <- seq(-6, 6, 2)
yseq <- seq(-2, 12, 2)
axis(1, xseq, labels=latex2exp::TeX(sprintf("$10^{%d}$", xseq)), 
     padj = -1,
     tck=-0.01,
     cex.axis=0.8, line=0.5)
axis(2, yseq, 
     las=1,
     tck=-0.01, hadj=0.5,
     labels=latex2exp::TeX(sprintf("$10^{%d}$", yseq)), 
     cex.axis=0.8, line=0.5)
legend(-6, 0, tools:::toTitleCase(levels(as.factor(d$group))), fill=group_cols,
       bty='n', border=0, cex=0.6, ncol=2)
title(xlab=latex2exp::TeX("mass (grams)"), line=2.2)
title(ylab=latex2exp::TeX("population density (individuals per $km^2$)"),   
      line=2.4)
dev.off()
par(opar)
