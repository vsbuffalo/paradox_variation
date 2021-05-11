# diversity_popsize_full.r -- 

library(tidyverse)
library(plotrix)
library(ggrepel)
library(mixtools)
library(scales)
library(sp)
library(grid)
library(latex2exp)

source('plot_funs.r')
source('color_scheme.r')

opar <- par(no.readonly=TRUE)

load('../../data/main_datasets.Rdata')
d <- da_dps
# summary(lm(pred_log10_N ~ log10_popsize, d))
# summary(lm(log10_diversity ~ log10_popsize, d))

# what to include
# x <- d$pred_log10_N  # the simplistic unregularized version
x <- d$log10_popsize
y <- d$diversity

set.seed(2)

xlims <- c(4, 16)

fit <- lm(y ~ x)

output <- TRUE
# output <- FALSE

if (output)
  pdf("diversity_popsize_linear.pdf", width=7, height=7/1.5)

par(mar=c(4, 4, 2, 2))

has_label <- d$label != ""
plot(x, y, col=phyla_cols[d$phylum], type='n', axes=FALSE,
     ylab='', xlab='', ylim=c(0, 0.09), xlim=xlims) 

# abline(a=log10(1e-8 * 4), b=1, lty=2, col='cornflowerblue')
# text(6, -0.9, latex2exp:::TeX("$\\mu = 10^{-8}$"), cex=0.5, col='cornflowerblue', xpd=TRUE)

# logN <- seq(3, 16, length.out=100)
# y1 <- log10(4e-8) + logN
# y2 <- log10(4e-10) + logN
# polygon(c(logN, rev(logN)), c(y1, rev(y2)), border=NA, col='gray88')

endp <- 11
logN <- seq(3, endp, length.out=100)
# y1 <- log10(4e-8) + logN
# y2 <- log10(4e-10) + logN
# polygon(c(logN, rev(logN)), c(y1, rev(y2)), border=NA, col='gray88')

# two allele
Kalleles <- 4
Kboth <- FALSE
palpha <- 0.4
if (Kalleles == 2 || Kboth) {
  y1 <- 10^(log10(4e-8) + logN - log10(2 + 2 * 4e-8 * 10^logN))
  y2 <- 10^(log10(4e-9) + logN - log10(2 + 2 * 4e-9 * 10^logN))
  text(6.2, 0.08, latex2exp:::TeX("$\\mu = 10^{-8}$"), 
       cex=0.5, col='gray52', xpd=TRUE)
  text(9.3, 0.08, latex2exp:::TeX("$\\mu = 10^{-9}$"), 
       cex=0.5, col='gray52', xpd=TRUE)
  polygon(c(logN, rev(logN)), c(y1, rev(y2)), border=NA, 
        col=alpha('gray88', palpha))
}
if (Kalleles == 4 || Kboth) {
  # four allele
  y1 <- 10^(log10(4e-8) + logN - log10(1 + 4/3 * 4e-8 * 10^logN))
  y2 <- 10^(log10(4e-9) + logN - log10(1 + 4/3 * 4e-9 * 10^logN))
  text(5.9, 0.09, latex2exp:::TeX("$\\mu = 10^{-8}$"), 
       cex=0.5, col='gray52', xpd=TRUE)
  text(8, 0.09, latex2exp:::TeX("$\\mu = 10^{-9}$"), 
       cex=0.5, col='gray52', xpd=TRUE)
  polygon(c(logN, rev(logN)), c(y1, rev(y2)), border=NA, 
        col=alpha('gray88', palpha))
}

# y <- predict(div_fit, newdata=data.frame(log10_popsize=logN))
# lines(logN, y, lty=2)
# abline(div_fit, lty=2, col='gray52')

points(x, y, pch=21, cex=1.3, lwd=0.4,
       bg=phyla_cols[d$phylum], col='white')
xseq <- seq(4, 16, 2)
axis(1, xseq, 
     labels=latex2exp::TeX(sprintf("$10^{%d}$", xseq)), 
     padj = -1,
     tck=-0.01,
     cex.axis=0.8, line=0.5)
yseq <- seq(0, 0.08, 0.02)
axis(2, yseq, 
     labels=yseq,
     las=1,
     tck=-0.01, hadj=0.5,
     cex.axis=0.8, line=0.5)

# rseq <- yseq - log10(4e-8)
# axis(4, yseq, 
#      labels=latex2exp::TeX(sprintf("$10^{%.2f}$", rseq)), 
#      las=1,
#      tck=-0.01, hadj=0.5,
#      cex.axis=0.8, line=0.5)
title(ylab=" pairwise diversity (differences per bp)", line=2.8, cex.lab=1.2)
title(xlab=latex2exp::TeX("approximate population size"), line=2.5, cex.lab=1.2)


# phyla <- as.factor(ld$phylum)
# for (i in 1:nphyla) {
#   l <- levels(phyla)[i]
#   z <- l == phyla
#   segments(min(ld$log10_size[z], na.rm=TRUE), -1.05 + 0.005*i, 
#            max(ld$log10_size[z], na.rm=TRUE), -1.05 + 0.005*i, 
#            col=phyla_cols[i], lwd=2)
# }
phyla_cols_s <- phyla_cols[sort(unique(d$phylum))]
legend(12.5, 0.075, names(phyla_cols_s), fill=phyla_cols_s,
       bty='n', border=0, cex=0.6, ncol=2)

if (output) dev.off()



par(opar)
