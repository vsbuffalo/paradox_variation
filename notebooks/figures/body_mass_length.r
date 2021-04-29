library(tidyverse)
library(scales)

opar <- par(no.readonly=TRUE)

pdf("body_mass_length.pdf", width=5, height=5)
par(mar=c(5, 5, 2, 1))

load('../../data/popsize_chains.Rdata')

pars <- rstan:::extract(bayes_fit, pars=c('alpha', 'beta'))
mean_intercept <- mean(pars$alpha)
mean_slope <- mean(pars$beta)
levels <- c("Annelida", "Arthropoda", "Chordata", 
            "Echinodermata", "Mollusca", "Nematoda", "NA")
d <- read_tsv('../../data/romiguier_et_al_2014_updated.tsv')   %>%
        mutate(log10_mass_g = log10(bodymass_g)) %>%
        filter_at(vars(log10_mass_g, log10_size), ~ !is.na(.)) %>%
        mutate(phylum = ifelse(is.na(phylum), "NA", phylum)) %>% 
        mutate(phylum = factor(phylum, levels=levels, ordered=TRUE))

ngroups <- length(unique(d$phylum)) - 1
group_cols <- c(hue_pal()(ngroups), 'gray52')

plot(d$log10_size, d$log10_mass_g, pch=19, 
     type='n', axes=FALSE, ylim=c(-4, 6),
     xlim=c(-3, 1),
     xlab='',
     ylab='')
points(d$log10_size, d$log10_mass_g, pch=21, 
       cex=1.5, lwd=0.4, bg=group_cols[as.factor(d$phylum)], col='white')
abline(a=mean_intercept, b=mean_slope, lty=2, col='gray52')

xseq <- seq(-3, 1, 1)
yseq <- seq(-4, 6, 2)
axis(1, xseq, labels=latex2exp::TeX(sprintf("$10^{%d}$", xseq)), 
     padj = -1,
     tck=-0.01,
     cex.axis=0.8, line=0.5)
axis(2, yseq, labels=latex2exp::TeX(sprintf("$10^{%d}$", yseq)), 
     las=1,
     tck=-0.01, hadj=0.5,
     cex.axis=0.8, line=0.5)
legend(-0.8, -2, tools:::toTitleCase(levels(as.factor(d$phylum))), fill=group_cols,
       bty='n', border=0, cex=0.6, ncol=2)
title(xlab=latex2exp::TeX("body length (meters)"), line=2.1)
title(ylab=latex2exp::TeX("body mass (grams)"), line=2.3)

dev.off()
par(opar)
