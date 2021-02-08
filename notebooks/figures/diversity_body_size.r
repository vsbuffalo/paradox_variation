# diversity_body_size.r -- create the body size / diversity figure

library(tidyverse)
library(plotrix)
library(ggrepel)
library(scales)
library(grid)
library(latex2exp)

source('plot_funs.r')

opar <- par(no.readonly=TRUE)

lds <- read_tsv('../../data/leffler_et_al_2012_updated.tsv') %>%
         mutate(log10_size=log10(size), 
                log10_diversity=log10(diversity))

set.seed(2)
p <- ggplot(lds) + 
  geom_point(aes(log10_size, log10_diversity, label=species, color=phylum), size=2) +
  geom_text_repel(aes(log10_size, log10_diversity, label=species, text=species), 
                  size=1.5, 
                  point.padding=0.29, 
                  #box.padding=0.2, 
                  max.iter=10000)

repel_pos <- extract_ggrepel(p)

# best to restart after messing with grid 
dev.off()


clean_species <- function(x) {
  cl <- sub('([^/]+)/.*', '\\1', x)
  sub('([^ ]+) ([^ ]+).*', '\\1 \\2', cl)
}

ld <- lds %>% filter(!is.na(log10_size), !is.na(log10_diversity)) %>%
         # simplify species names, remove the species1/species2 in Leffler
         mutate(species = map_chr(species, clean_species))


nphyla <- length(unique(ld$phylum))
phyla_cols <- hue_pal()(nphyla)


pdf("diversity_body_size.pdf", width=7.3, height=7.3/1.5)
par(mar=c(4, 4, 2, 1))
x <- ld$log10_size
y <- ld$log10_diversity
tx <- repel_pos[, 1]
ty <- repel_pos[, 2]
plot(x, y, col=phyla_cols[as.factor(ld$phylum)], type='n', axes=FALSE,
     ylab='', xlab='', xlim=c(-5.3, 2))
adjust <- 0
segments(x, y, tx, ty + adjust*ifelse(ty > y, -0.02, 0.02), col=alpha('gray42', 0.5), lwd=0.5)
points(x, y, pch=21, cex=1.3, lwd=0.4,
       bg=phyla_cols[as.factor(ld$phylum)], col='gray32')
#text(tx, ty, ld$species, cex=0.26, border=FALSE)
boxed.labels(tx, ty, ld$species, cex=0.29, border=FALSE, bg=alpha('white', 0.6), xpad=1, ypad=1)
axis(1, seq(-5, 2), labels=latex2exp::TeX(sprintf("$10^{%d}$", seq(-5, 2))), 
     cex.axis=0.8, line=0.5)
axis(2, seq(-4, -1), labels=latex2exp::TeX(sprintf("$10^{%d}$", seq(-4, -1))),
     cex.axis=0.8, line=0.5)
title(ylab="pairwise diveristy", line=2.8, cex.lab=1.2)
title(xlab="body size (meters)", line=2.8, cex.lab=1.2)


# phyla <- as.factor(ld$phylum)
# for (i in 1:nphyla) {
#   l <- levels(phyla)[i]
#   z <- l == phyla
#   segments(min(ld$log10_size[z], na.rm=TRUE), -1.05 + 0.005*i, 
#            max(ld$log10_size[z], na.rm=TRUE), -1.05 + 0.005*i, 
#            col=phyla_cols[i], lwd=2)
# }
legend(-5.4, -3.6, levels(as.factor(ld$phylum)), fill=phyla_cols,
       bty='n', border=0, cex=0.6, ncol=3)
dev.off()

par(opar)

