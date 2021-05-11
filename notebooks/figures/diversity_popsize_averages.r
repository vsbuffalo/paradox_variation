# diversity_popsize_averages.r -- 

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
load('../../data/phyla_diversity_popsize_ols.Rdata')

# alias this -- less typing
d <- da_dps

# what to include
# x <- d$pred_log10_N  # the simplistic unregularized version
x <- d$log10_popsize
y <- d$log10_diversity

output <- TRUE

if (output)
  pdf("diversity_popsize_averages.pdf", width=7, height=7/1.5)

par(mar=c(4, 4, 2, 2))

xlims <- c(4, 16)
has_label <- d$label != ""
plot(x, y, col=phyla_cols[d$phylum], type='n', axes=FALSE,
     ylab='', xlab='', ylim=c(-4.1, -0.9), xlim=xlims) 

# text(6, -0.9, latex2exp:::TeX("$\\mu = 10^{-8}$"), cex=0.5, col='gray52', xpd=TRUE)
# text(9, -0.9, latex2exp:::TeX("$\\mu = 10^{-10}$"), cex=0.5, col='gray52', xpd=TRUE)

logN <- seq(3, 16, length.out=100)
y1 <- log10(4e-8) + logN
y2 <- log10(4e-10) + logN
# polygon(c(logN, rev(logN)), c(y1, rev(y2)), border=NA, col='gray88')

fits <- ols_phyla$ols
names(fits) <- ols_phyla$phylum

for (i in seq_along(fits)) {
  lps <- ols_phyla$data[[i]]$log10_popsize
  zz <- seq(min(lps), max(lps), length.out=100)
  yy <- coef(fits[[i]])[1] + coef(fits[[i]])[2]*zz
  lines(zz, yy, lwd=1.2, lty=2, 
         col=phyla_cols[names(fits)[i]])
}

adjust <- 0
xc <- x[has_label]
yc <- y[has_label]
points(x, y, pch=21, cex=1.3, lwd=0.4,
       bg=phyla_cols[d$phylum], col='white')
xseq <- seq(4, 16, 2)
axis(1, xseq, 
     labels=latex2exp::TeX(sprintf("$10^{%d}$", xseq)), 
     padj = -1,
     tck=-0.01,
     cex.axis=0.8, line=0.5)
yseq <- seq(-4, -1)
axis(2, yseq, 
     labels=latex2exp::TeX(sprintf("$10^{%d}$", yseq)), 
     las=1,
     tck=-0.01, hadj=0.5,
     cex.axis=0.8, line=0.5)

# rseq <- yseq - log10(4e-8)
# axis(4, yseq, 
#      labels=latex2exp::TeX(sprintf("$10^{%.2f}$", rseq)), 
#      las=1,
#      tck=-0.01, hadj=0.5,
#      cex.axis=0.8, line=0.5)
title(ylab="   pairwise diversity (differences per bp)", 
      line=2.8, cex.lab=1.2)
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
legend(10, -3.6, names(phyla_cols_s), fill=phyla_cols_s,
       bty='n', border=0, cex=0.6, ncol=3)

if (output) dev.off()



par(opar)
