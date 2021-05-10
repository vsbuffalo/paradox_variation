# species_counts.r

library(tidyverse)
source('color_scheme.r')

d <- read_tsv('../../data/species_by_class.tsv') %>% 
  filter(is.finite(frac))
opar <- par(no.readonly=TRUE)

size <- 7
pdf('species_counts.pdf', width=size*0.773, size)
par(mar=c(5, 0.1, 0.1, 0.6))
x <- seq_along(d$class)
plot(log10(d$frac), x, type='n', axes=FALSE, ylab='', xlab='', 
     xlim=c(-6.5, -0.5))
axis(1, seq(-1, -5), latex2exp::TeX(sprintf("$10^{%d}$", seq(-1, -5))))

cols <- match(d$phylum, names(phyla_cols))
points(log10(d$frac), x, pty='n', cex=1.5*(d$nspecies)^0.35, pch=21,
       col=phyla_cols[cols], lwd=3, bg=scales::alpha(phyla_cols[cols], 0.6))
text(log10(d$frac) + 0.11*d$nspecies^0.3, x, d$class, 
     xpd=TRUE, srt=0, cex=0.8, adj=c(0, NA))

text(log10(d$frac), x, d$nspecies, cex=0.6, adj=c(0.5, NA))

z = factor(d$phylum)
phy <- as.integer(z)
for (i in 1:max(phy)) {
  idx <- which(phy == i) 
  text(-5.5, mean(idx), levels(z)[i], cex=0.8, adj=c(1, NA))
  segments(-5.4, min(idx) - 0.1, -5.4, max(idx) + 0.1, lwd=3, lend=2, 
          col=phyla_cols[i])
}

title(xlab='                   fraction of estimated total species included')

dev.off()

# phyla_cols_s <- phyla_cols[unique(d$phylum)]
# legend(-5, 15, names(phyla_cols_s), fill=phyla_cols_s,
#        bty='n', border=0, cex=0.7, ncol=1)
par(opar)
