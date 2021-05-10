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

# for stan fit
load('../../data/diversity_popsize_chains.Rdata')
load('../../data/main_datasets.Rdata')
phi <- mean(rstan:::extract(dps_fit, pars='phi')$phi)
beta <- mean(rstan:::extract(dps_fit, pars='beta')$beta)

# alias this -- less typing
d <- da_dps %>% filter(species != 'Equus ferus przewalskii')

#d = read_tsv('../../data/combined_data.tsv')
#d %>% filter(class == 'Aves') %>% select(species, log10_range, log10_popsize)
 
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
keep_sps <- c("Homo sapiens", "Apis mellifera", "Canis lupus", "Caenorhabditis elegans")
drosoph <- c(
             "Drosophila americana", 
             "Drosophila ananassae", 
             #"Drosophila arizonae",
             #"Drosophila bipectinata", 
             #"Drosophila buzzatii",  
             #"Drosophila kikkawai",
             #"Drosophila malerkotliana", 
             "Drosophila melanogaster", 
             "Drosophila miranda",
             "Drosophila mojavensis", 
             #"Drosophila novamexicana", 
             #"Drosophila parabipectinata",
             "Drosophila persimilis", 
             "Drosophila pseudoobscura", 
             "Drosophila sechellia",
             "Drosophila simulans", 
             "Drosophila subobscura")

label <- (label | label %in% drosoph) | label %in% keep_sps
# contained <- rep(TRUE, length(x))

d$label <- ifelse(label, d$species, "")

# ignore a cluster of points a little over Nc=10^8
# which has human
#ignore <- c('Aquila clanga', 'Reticulitermes lucifugus', 'Reticulitermes grassei')
#d$label <- ifelse(d$label %in% ignore, "", d$label)


# simplify some species name (Leffler data has combined them in some cases)
d$label[d$species == 'Gorilla beringei/gorilla'] <- "Gorilla gorilla"
d$label[d$species == 'Pongo abelii/pygmaeus'] <- "Pongo abelii"
d$label[d$species == 'Nomascus gabriellae/leucogenys'] <- "Nomascus gabriellae"
d$label[d$species == 'Aquila clanga/pomarina'] <- "Aquila clanga"

set.seed(2)

xlims <- c(4, 16)
ylims <- c(-4.1, 0)


REPEL_CACHED_FILE <- 'repel_full.tsv'
FORCE <- TRUE
#FORCE <- FALSE
if (FORCE || !file.exists(REPEL_CACHED_FILE)) {
  p <- ggplot(d) + 
    geom_point(aes(log10_popsize, log10_diversity, 
                   label=species, color=phylum), size=2) +
    geom_text_repel(data=d,
                    mapping=aes(log10_popsize, log10_diversity, 
                                label=label, text=species), 
                                size=2.1,
                                point.padding=10, 
                                force = 1.5,
                                min.segment.length=0.5,
                                seed=1,
                                box.padding=1e-1, 
                                ylims=ylims,
                                max.iter=1000000, xlim=xlims) +
    xlim(xlims[1], xlims[2]) + theme(legend.position = "none") + 
    ylim(ylims[1], ylims[2])
  
  repel_pos <- extract_ggrepel(p)
  write_tsv(repel_pos, path=REPEL_CACHED_FILE)
} else {
  repel_pos <- read_tsv(REPEL_CACHED_FILE)
}

# best to restart after messing with grid 
dev.off()


fit <- lm(y ~ x)

output <- TRUE
# output <- FALSE

if (output)
  pdf("diversity_popsize_full.pdf", width=7, height=7/1.5)

par(mar=c(4, 4, 1, 2))

tx <- repel_pos$x
ty <- repel_pos$y
tl <- repel_pos$label

td <- repel_pos %>% 
       left_join(d %>% select(label, x=log10_popsize, y=log10_diversity),
                 by='label', suffix=c('_text', '_point')) %>%
       as_tibble() %>% filter(is.finite(x_point), is.finite(y_point))

print(sprintf("number of labels: %d", nrow(td)))

plot(x, y, col=phyla_cols[d$phylum], type='n', axes=FALSE,
     ylab='', xlab='', ylim=ylims, xlim=xlims) 

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
  y1 <- log10(4e-8) + logN - log10(2 + 2 * 4e-8 * 10^logN)
  y2 <- log10(4e-9) + logN - log10(2 + 2 * 4e-9 * 10^logN)
  lines(xx <- seq(endp-0.1, 18), rep(log10(1/2), length(xx)), 
        col=alpha('gray88', palpha), lty=2, lwd=1, lend=1)

  text(7.2, -0.5, latex2exp:::TeX("$\\mu = 10^{-8}$"), 
       cex=0.5, col='gray52', xpd=TRUE)
  text(10.2, -0.5, latex2exp:::TeX("$\\mu = 10^{-9}$"), 
       cex=0.5, col='gray52', xpd=TRUE)
  polygon(c(logN, rev(logN)), c(y1, rev(y2)), border=NA, 
        col=alpha('gray88', palpha))
}
if (Kalleles == 4 || Kboth) {
  # four allele
  y1 <- log10(4e-8) + logN - log10(1 + 4/3 * 4e-8 * 10^logN)
  y2 <- log10(4e-9) + logN - log10(1 + 4/3 * 4e-9 * 10^logN)
  lines(xx <- seq(endp-0.1, 18), rep(log10(3/4), length(xx)), 
        col=alpha('gray88', 0.4), lty=2, lwd=1, lend=1)
  text(7.1, -0.3, latex2exp:::TeX("$\\mu = 10^{-8}$"), 
       cex=0.5, col='gray52', xpd=TRUE)
  text(9.2, -0.3, latex2exp:::TeX("$\\mu = 10^{-9}$"), 
       cex=0.5, col='gray52', xpd=TRUE)
  polygon(c(logN, rev(logN)), c(y1, rev(y2)), border=NA, 
        col=alpha('gray88', palpha))
}

# text(6.7, -0.3, latex2exp:::TeX("$\\mu = 10^{-8}$"), cex=0.5, col='gray52', xpd=TRUE)
# text(10.5, -0.3, latex2exp:::TeX("$\\mu = 10^{-10}$"), cex=0.5, col='gray52', xpd=TRUE)


# y <- predict(div_fit, newdata=data.frame(log10_popsize=logN))
# lines(logN, y, lty=2)
# abline(dps_ols, lty=2, col='gray52')
z <- seq(xlims[1], xlims[2], length.out=100)
lines(z, phi + beta*z, lty=2, col='cornflowerblue')
lines(z, coef(dps_ols)[1] + z*coef(dps_ols)[2], lty=2, col='gray52')

adjust <- 0
with(td,
     segments(x_point, y_point, x_text, y_text, col=alpha('gray42', 0.5), lwd=0.5))
points(x, y, pch=21, cex=1.3, lwd=0.4,
       bg=phyla_cols[d$phylum], col='white')
with(td, text(x_text, y_text, label, cex=0.29))
# boxed.labels(tx, ty, d$label[has_label], cex=0.5, border=FALSE, xpad=1, ypad=1)
xseq <- seq(4, 16, 2)
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

# rseq <- yseq - log10(4e-8)
# axis(4, yseq, 
#      labels=latex2exp::TeX(sprintf("$10^{%.2f}$", rseq)), 
#      las=1,
#      tck=-0.01, hadj=0.5,
#      cex.axis=0.8, line=0.5)
title(ylab="   pairwise diversity (differences per bp)", line=2.8, cex.lab=1.2)
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
