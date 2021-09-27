library(tidyverse)

opar <- par(no.readonly=TRUE)
drl <- read_tsv('../../data/redlist_popsizes.tsv') %>%
        mutate(log10_range = log10(range), log10_eoo_km2=log10(eoo_km2))

pdf('iucn_redlist_ranges.pdf', width=5, height=5)

fit <- lm(log10_range ~ log10_eoo_km2, data=drl)
par(mar=c(5, 5, 2, 2))
plot(log10_range ~ log10_eoo_km2, data=drl, pch=19,
     axes=FALSE, ylab='',
     xlab='')
xlabs <- seq(4, 8)
abline(b=1, a=0, col='gray42', lty=2)
abline(fit)
axis(1, xlabs, labels=latex2exp::TeX(sprintf("$10^{%d}$", xlabs)),
     padj = -1,
     tck=-0.01,
     cex.axis=0.8, line=0.5)
axis(2, xlabs, labels=latex2exp::TeX(sprintf("$10^{%d}$", xlabs)),
     las=1,
     tck=-0.01, hadj=0.5,
     cex.axis=0.8, line=0.5)
corr <- cor(drl$log10_eoo_km2, drl$log10_range, use='pairwise.complete')
corpval <- cor.test(drl$log10_eoo_km2, drl$log10_range, use='pairwise.complete')$p.value
xpos <- 6.2
ypos <- 4.5
text(xpos, ypos + .2, latex2exp::TeX(sprintf('$\\rho = %.2f$', round(corr, 2))), 
     cex=0.8, adj=0)

mag <- round(log10(corpval)) - 1
mantissa <- corpval * 10^abs(mag)
text(xpos, ypos, latex2exp::TeX(sprintf('p-value = $%.2f \\times 10^{%d}$', 
                                  mantissa, mag)), cex=0.8, adj=0)

title(ylab=latex2exp:::TeX('estimated range ($km^2$)'), line=2.8) 
title(xlab=latex2exp:::TeX('IUCN Red List estimated extent of occurrence ($km^2$)'),
      line=2.4) 
dev.off()
par(opar)
