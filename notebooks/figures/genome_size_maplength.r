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
      mutate(log10_genome_size = log10(genome_size),
             log10_map_length = log10(map_length)) %>% 
      filter_at(vars(log10_genome_size, map_length), 
                ~ is.finite(.)) %>% 
      filter(species != "Ambystoma tigrinum")


# what to include
# x <- d$pred_log10_N  # the simplistic unregularized version
x <- d$log10_genome_size
y <- d$log10_map_length

plot(x, y)

X <- cbind(x, y)
center <- colMeans(X)
sigma <- cov(X)
coords <- ellipse(mu=center, sigma=sigma, alpha=0.4, draw=FALSE)

label <- !as.logical(point.in.polygon(x, y, coords[, 1], coords[, 2]))
# contained <- rep(TRUE, length(x))

d$label <- ifelse(label, d$species, "")

set.seed(2)

xlims <- c(-1.2, 1)

REPEL_CACHED_FILE <- 'repel_genome_size_maplength.tsv'
FORCE <- TRUE
if (FORCE || !file.exists(REPEL_CACHED_FILE)) {
  p <- ggplot(d) + 
    geom_point(aes(log10_genome_size, log10_map_length, 
                   color=phylum), size=2) +
    geom_text_repel(data=d,
                    mapping=aes(log10_genome_size, 
                                log10_map_length, 
                                label=label),
                                size=2.1, 
                                point.padding=10, 
                                force=2,
                                min.segment.length=0.1,
                                seed=1,
                                box.padding=1.3e-1, 
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
  pdf("genome_size_maplength.pdf", width=7, height=7/1.5)

par(mar=c(4, 4, 2, 2))

tx <- repel_pos$x
ty <- repel_pos$y
tl <- repel_pos$label

td <- repel_pos %>% 
       left_join(d %>% select(label, y=map_length, x=log10_genome_size),
      by='label', suffix=c('_text', '_point')) %>%
       as_tibble() %>% filter(is.finite(x_point), is.finite(y_point))

plot(x, y, col=phyla_cols[d$phylum], 
     type='n', axes=FALSE, ylim=c(-0.3, 2),
     ylab='', xlab='', xlim=xlims) 

logN <- seq(3, 16, length.out=100)


# y <- predict(div_fit, newdata=data.frame(log10_popsize=logN))
# lines(logN, y, lty=2)
abline(fit, lty=2, lwd=1.2, col=alpha('gray52', 0.5))

dfit <- d %>% group_by(phylum) %>%
     nest() %>%
     mutate(n = map_int(data, nrow)) %>%
     filter(n > 5) %>%
     mutate(fit = map(data, ~ lm(log10_map_length ~ log10_genome_size, .)))

ave_lines <- TRUE 
if (ave_lines) {
  for (i in seq_along(dfit$fit)) { 
    fit <- dfit$fit[[i]]
    phylum <- dfit$phylum[[i]]
    xs <- sort(dfit$data[[i]]$log10_genome_size)
    ys <- predict(fit, newdata=data.frame(log10_genome_size=xs))
    lines(xs, ys, lty=2, lwd=1.2, col=phyla_cols[phylum])
  }
}


with(td,
     segments(x_point, y_point, x_text, y_text, col=alpha('gray42', 0.5), lwd=0.5))
points(x, y, pch=21, cex=1.3, lwd=0.4,
       bg=phyla_cols[d$phylum], col='white')
with(td, text(x_text, y_text, label, cex=0.29))
# boxed.labels(tx, ty, d$label[has_label], cex=0.5, border=FALSE, xpad=1, ypad=1)
xseq <- seq(-1, 1, 1)
axis(1, xseq, 
     labels=latex2exp::TeX(sprintf("$10^{%d}$", xseq)), 
     padj = -1,
     tck=-0.01,
     cex.axis=0.8, line=0.5)
yseq <- seq(0, 2, 1)
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
title(xlab="genome size (Gb)", line=2.8, cex.lab=1.2)
title(ylab=latex2exp::TeX("map length (Morgans)"), line=2.5, cex.lab=1.2)


# phyla <- as.factor(ld$phylum)
# for (i in 1:nphyla) {
#   l <- levels(phyla)[i]
#   z <- l == phyla
#   segments(min(ld$log10_size[z], na.rm=TRUE), -1.05 + 0.005*i, 
#            max(ld$log10_size[z], na.rm=TRUE), -1.05 + 0.005*i, 
#            col=phyla_cols[i], lwd=2)
# }
phyla_cols_s <- phyla_cols[unique(d$phylum)]
legend(-0.1, 0, names(phyla_cols_s), fill=phyla_cols_s,
       bty='n', border=0, cex=0.6, ncol=3)

if (output) dev.off()



par(opar)
