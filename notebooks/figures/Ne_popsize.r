# Ne_popsize.r -- 

library(tidyverse)
library(plotrix)
library(ggrepel)
library(mixtools)
library(scales)
library(sp)
library(grid)
library(latex2exp)

source('plot_funs.r')

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
y <- d$log10_diversity - log10(4e-8)

plot(x, y)

X <- cbind(x, y)
center <- colMeans(X)
sigma <- cov(X)

coords <- ellipse(mu=center, sigma=sigma, alpha=0.4, draw=FALSE)

label <- !as.logical(point.in.polygon(x, y, coords[, 1], coords[, 2]))
# contained <- rep(TRUE, length(x))

d$label <- ifelse(label, d$species, "")

set.seed(2)

xlims <- c(4, 16.5)

REPEL_CACHED_FILE <- 'Ne_repel.tsv'
if (!file.exists(REPEL_CACHED_FILE)) {
  p <- ggplot(d) + 
    geom_point(aes(log10_popsize, log10_diversity - log10(4e-8), 
                   label=species, color=phylum), size=2) +
    geom_text_repel(data=d,
                    mapping=aes(log10_popsize, log10_diversity - log10(4e-8), 
                                label=label, text=species), 
                                size=1.5, 
                                point.padding=0.6, 
                                force = 2,
                                min.segment.length=0.4,
                                seed=1,
                                box.padding=0.3, 
                                max.iter=100000, xlim=xlims) +
    xlim(xlims[1], xlims[2]) + theme(legend.position = "none")
  
  repel_pos <- extract_ggrepel(p)
  write_tsv(repel_pos, path=REPEL_CACHED_FILE)
} else {
  repel_pos <- read_tsv(REPEL_CACHED_FILE)
}

# best to restart after messing with grid 
dev.off()

nphyla <- length(unique(d$phylum))
phyla_cols <- hue_pal()(nphyla)

fit <- lm(y ~ x)

output <- TRUE

if (output)
  pdf("Ne_popsize.pdf", width=7, height=7/1.5)

par(mar=c(4, 4, 2, 4))

has_label <- d$label != ""
tx <- repel_pos$x
ty <- repel_pos$y
ymax <- 6.3
plot(x, y, col=phyla_cols[as.factor(d$phylum)], type='n', axes=FALSE,
     ylab='', xlab='', ylim=c(3, ymax), xlim=xlims) 


# abline(a=log10(1e-8 * 4), b=1, lty=2, col='cornflowerblue')
# text(6, -0.9, latex2exp:::TeX("$\\mu = 10^{-8}$"), cex=0.5, col='cornflowerblue', xpd=TRUE)
# abline(a=log10(1e-10 * 4), b=1, lty=2, col='gray52')

logN <- seq(3, 16, length.out=100)
logN2 <- seq(3, 8.35, length.out=100)
logN3 <- seq(3, 6.35, length.out=100)

# y <- predict(div_fit, newdata=data.frame(log10_popsize=logN))
# lines(logN, y, lty=2)
# abline(div_fit, lty=2, col='gray52')

lines(logN3, logN3, col='gray52', cex=0.5, lty=2)
lines(logN2, -log10(4e-8) + log10(4e-10) + logN2, col='gray52', cex=0.5, lty=2)
text(4.1, 3.6, latex2exp:::TeX("$\\mu = 10^{-8}$"), cex=0.5, col='gray52', xpd=TRUE)
text(6.1, 3.6, latex2exp:::TeX("$\\mu = 10^{-10}$"), cex=0.5, col='gray52', xpd=TRUE)
# abline(a=0, b=1)
# abline(a=-log10(4e-8) + log10(4e-10), b=1)

adjust <- 0
xc <- x[has_label]
yc <- y[has_label]
segments(xc, yc, tx, ty, col=alpha('gray42', 0.5), lwd=0.5)
points(x, y, pch=21, cex=1.3, lwd=0.4,
       bg=phyla_cols[as.factor(d$phylum)], col='white')
text(tx, ty, d$label[has_label], cex=0.3)
# boxed.labels(tx, ty, d$label[has_label], cex=0.5, border=FALSE, xpad=1, ypad=1)
xseq <- seq(4, 16, 2)
axis(1, xseq, 
     labels=latex2exp::TeX(sprintf("$10^{%d}$", xseq)), 
     padj = -1,
     tck=-0.01,
     cex.axis=0.8, line=0.5)
# axis(2, yseq, 
#      labels=latex2exp::TeX(sprintf("$10^{%d}$", yseq)), 
#      las=1,
#      tck=-0.01, hadj=0.5,
#      cex.axis=0.8, line=0.5)

yseq <- seq(3, floor(ymax))
axis(2, yseq, 
     labels=latex2exp::TeX(sprintf("$10^{%d}$", yseq)), 
     las=1,
     tck=-0.01, hadj=0.5,
     cex.axis=0.8, line=0.5)

yseq2 <- yseq + log10(4e-8 / 4e-10)
axis(4, yseq, 
     labels=latex2exp::TeX(sprintf("$10^{%d}$", yseq2)), 
     las=1,
     tck=-0.01, hadj=0.5,
     cex.axis=0.8, line=0.5)



# title(ylab="pairwise diversity", line=2.8, cex.lab=1.2)
title(xlab=latex2exp::TeX("approximate population size"), line=2.5, cex.lab=1.3)
mtext(latex2exp::TeX("approximate $N_e$ assuming $\\mu = 10^{-8}$"), 2, line=2.5, cex.lab=1.2)
text(18.5, 4.5, latex2exp::TeX("approximate $N_e$ assuming $\\mu = 10^{-10}$"), 
     srt=-90, xpd=TRUE)


# phyla <- as.factor(ld$phylum)
# for (i in 1:nphyla) {
#   l <- levels(phyla)[i]
#   z <- l == phyla
#   segments(min(ld$log10_size[z], na.rm=TRUE), -1.05 + 0.005*i, 
#            max(ld$log10_size[z], na.rm=TRUE), -1.05 + 0.005*i, 
#            col=phyla_cols[i], lwd=2)
# }
legend(9.5, 3.5, levels(as.factor(d$phylum)), fill=phyla_cols,
       bty='n', border=0, cex=0.6, ncol=3)

if (output) dev.off()



par(opar)
