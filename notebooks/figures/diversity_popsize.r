# diversity_popsize.r -- 

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

# for div_fit
load('../../data/lognorm_fulldata_with_social.Rdata')


d <- read_tsv('../../data/combined_data.tsv') %>%
      filter(kingdom == 'Animalia') %>%
      mutate(log10_diversity = log10(diversity)) %>% 
      filter_at(vars(log10_diversity, pred_log10_N, log10_popsize), 
                ~ is.finite(.))


# summary(lm(pred_log10_N ~ log10_popsize, d))
# summary(lm(log10_diversity ~ log10_popsize, d))

# what to include
# x <- d$pred_log10_N  # the simplistic unregularized version
x <- d$log10_popsize
y <- d$log10_diversity

plot(x, y)

X <- cbind(x, y)
center <- colMeans(X)
sigma <- cov(X)

coords <- ellipse(mu=center, sigma=sigma, alpha=0.35, draw=FALSE)

label <- !as.logical(point.in.polygon(x, y, coords[, 1], coords[, 2]))
# add in all drosophila
keep_sps <- c("Homo sapiens", "Apis mellifera", "Canis lupus")
label <- (label | d$family == 'Drosophilidae') | label %in% keep_sps
# contained <- rep(TRUE, length(x))

d$label <- ifelse(label, d$species, "")

set.seed(2)

xlims <- c(4, 16.5)

REPEL_CACHED_FILE <- 'repel.tsv'
if (!file.exists(REPEL_CACHED_FILE)) {
  p <- ggplot(d) + 
    geom_point(aes(log10_popsize, log10_diversity, 
                   label=species, color=phylum), size=2) +
    geom_text_repel(data=d,
                    mapping=aes(log10_popsize, log10_diversity, 
                                label=label, text=species), 
                                size=1.5, 
                                point.padding=0.45, 
                                force = 2,
                                min.segment.length=0.3,
                                seed=1,
                                box.padding=0.2, 
                                max.iter=1000000, xlim=xlims) +
    xlim(xlims[1], xlims[2]) + theme(legend.position = "none")
  
  repel_pos <- extract_ggrepel(p)
  write_tsv(repel_pos, path=REPEL_CACHED_FILE)
} else {
  repel_pos <- read_tsv(REPEL_CACHED_FILE)
}

# best to restart after messing with grid 
dev.off()


fit <- lm(y ~ x)

output <- TRUE

if (output)
  pdf("diversity_popsize.pdf", width=7, height=7/1.5)

par(mar=c(4, 4, 2, 2))

has_label <- d$label != ""
tx <- repel_pos$x
ty <- repel_pos$y
plot(x, y, col=phyla_cols[d$phylum], type='n', axes=FALSE,
     ylab='', xlab='', ylim=c(-4.1, -0.9), xlim=xlims) 

# abline(a=log10(1e-8 * 4), b=1, lty=2, col='cornflowerblue')
# text(6, -0.9, latex2exp:::TeX("$\\mu = 10^{-8}$"), cex=0.5, col='cornflowerblue', xpd=TRUE)
text(6.1, -0.9, latex2exp:::TeX("$\\mu = 10^{-8}$"), cex=0.5, col='gray52', xpd=TRUE)
# abline(a=log10(1e-10 * 4), b=1, lty=2, col='gray52')
text(9.1, -0.9, latex2exp:::TeX("$\\mu = 10^{-10}$"), cex=0.5, col='gray52', xpd=TRUE)

# logN <- seq(3, 16, length.out=100)
# y1 <- log10(4e-8) + logN
# y2 <- log10(4e-10) + logN
# polygon(c(logN, rev(logN)), c(y1, rev(y2)), border=NA, col='gray88')

logN <- seq(3, 16, length.out=100)
# y1 <- log10(4e-8) + logN
# y2 <- log10(4e-10) + logN
# polygon(c(logN, rev(logN)), c(y1, rev(y2)), border=NA, col='gray88')

y1 <- log10(4e-8) + logN - log10(1 + 4/3 * 4e-8 * 10^logN)
y2 <- log10(4e-10) + logN - log10(1 + 4/3 * 4e-10 * 10^logN)
polygon(c(logN, rev(logN)), c(y1, rev(y2)), border=NA, col='gray88')
lines(xx <- seq(11.5, 17), rep(log10(3/4), length(xx)), col='gray88', lty=2, lwd=1, lend=1)

# text(6.7, -0.3, latex2exp:::TeX("$\\mu = 10^{-8}$"), cex=0.5, col='gray52', xpd=TRUE)
# text(10.5, -0.3, latex2exp:::TeX("$\\mu = 10^{-10}$"), cex=0.5, col='gray52', xpd=TRUE)


# y <- predict(div_fit, newdata=data.frame(log10_popsize=logN))
# lines(logN, y, lty=2)
abline(div_fit, lty=2, col='gray52')

adjust <- 0
xc <- x[has_label]
yc <- y[has_label]
segments(xc, yc, tx, ty, col=alpha('gray42', 0.5), lwd=0.5)
points(x, y, pch=21, cex=1.3, lwd=0.4,
       bg=phyla_cols[d$phylum], col='white')
text(tx, ty, d$label[has_label], cex=0.3)
# boxed.labels(tx, ty, d$label[has_label], cex=0.5, border=FALSE, xpad=1, ypad=1)
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
title(ylab="pairwise diversity", line=2.8, cex.lab=1.2)
title(xlab=latex2exp::TeX("approximate population size"), line=2.5, cex.lab=1.2)


# phyla <- as.factor(ld$phylum)
# for (i in 1:nphyla) {
#   l <- levels(phyla)[i]
#   z <- l == phyla
#   segments(min(ld$log10_size[z], na.rm=TRUE), -1.05 + 0.005*i, 
#            max(ld$log10_size[z], na.rm=TRUE), -1.05 + 0.005*i, 
#            col=phyla_cols[i], lwd=2)
# }
phyla_cols_s <- phyla_cols[unique(d$phylum)]
legend(10, -3.6, names(phyla_cols_s), fill=phyla_cols_s,
       bty='n', border=0, cex=0.6, ncol=3)

if (output) dev.off()



par(opar)
