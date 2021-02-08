## NOTE:
# if for some reason repel isn't working, the xlim is likely
# messedup
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

d <- read_tsv('../../data/combined_data.tsv') %>%
      mutate(log10_genome_size = log10(genome_size)) %>% 
      filter_at(vars(log10_genome_size, log10_popsize), 
                ~ is.finite(.)) %>% 
      filter(species != "Ambystoma tigrinum")


# what to include
# x <- d$pred_log10_N  # the simplistic unregularized version
x <- d$log10_popsize
y <- d$log10_genome_size

plot(x, y)

X <- cbind(x, y)
center <- colMeans(X)
sigma <- cov(X)
coords <- ellipse(mu=center, sigma=sigma, alpha=0.3,
                  draw=FALSE)

label <- !as.logical(point.in.polygon(x, y, coords[, 1], coords[, 2]))
# contained <- rep(TRUE, length(x))

d$label <- ifelse(label, d$species, "")

set.seed(2)

xlims <- c(3.9, 18)

REPEL_CACHED_FILE <- 'repel_genome_size_popsize.tsv'
# FORCE <- TRUE
FORCE <- FALSE
if (FORCE || !file.exists(REPEL_CACHED_FILE)) {
  p <- ggplot(d) + 
    geom_point(aes(log10_popsize, log10_genome_size, 
                   label=species, color=phylum), size=2) +
    geom_text_repel(data=d,
                    mapping=aes(log10_popsize, log10_genome_size, 
                                label=label, text=species),
                                size=1.5, 
                                point.padding=0.6, 
                                force=2,
                                min.segment.length=0.4,
                                seed=1,
                                box.padding=0.6, 
                                max.iter=100000, xlim=xlims) +
    xlim(xlims[1], xlims[2]) + theme(legend.position = "none")
  repel_pos <- extract_ggrepel(p)
  write_tsv(repel_pos, path=REPEL_CACHED_FILE)
} else {
  repel_pos <- read_tsv(REPEL_CACHED_FILE)
}

# best to restart after messing with grid 
# dev.off()

fit <- lm(y ~ x)

output <- TRUE

if (output)
  pdf("genome_size_popsize.pdf", width=7, height=7/1.5)

par(mar=c(4, 4, 2, 2))

has_label <- d$label != ""
tx <- repel_pos$x
ty <- repel_pos$y
plot(x, y, col=phyla_cols[d$phylum], 
     type='n', axes=FALSE, ylim=c(-1.2, 1),
     ylab='', xlab='', xlim=xlims) 

logN <- seq(3, 16, length.out=100)


# y <- predict(div_fit, newdata=data.frame(log10_popsize=logN))
# lines(logN, y, lty=2)
abline(fit, lty=2, lwd=1.2, col=alpha('gray52', 0.5))

dfit <- d %>% group_by(phylum) %>%
     nest() %>%
     mutate(n = map_int(data, nrow)) %>%
     filter(n > 5) %>%
     mutate(fit = map(data, ~ lm(log10_genome_size ~ log10_popsize, .)))

ave_lines <- FALSE
if (ave_lines) {
  for (i in seq_along(dfit$fit)) { 
    fit <- dfit$fit[[i]]
    phylum <- dfit$phylum[[i]]
    xs <- sort(dfit$data[[i]]$log10_popsize)
    ys <- predict(fit, newdata=data.frame(log10_popsize=xs))
    lines(xs, ys, lty=2, lwd=1.2, col=phyla_cols[phylum])
  }
}

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
yseq <- seq(-1, 1)
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
title(ylab="genome size (Gb)", line=2.8, cex.lab=1.2)
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
legend(4, -1, names(phyla_cols_s), fill=phyla_cols_s,
       bty='n', border=0, cex=0.6, ncol=3)

if (output) dev.off()



par(opar)
