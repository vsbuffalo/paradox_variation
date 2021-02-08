# map_length_body_size.r -- create the body size / map length figure

library(tidyverse)
library(ggrepel)
library(viridis)
library(wesanderson)
library(shape)


load('../../data/body_size_map_length.Rdata')
opar <- par(no.readonly=TRUE)
arrow_col <- 'gray32'

pdf <- TRUE
# pdf <- FALSE
if (pdf) pdf("map_length_body_size_nls.pdf", width=7.3, height=7.3)


#par(mfrow=c(1, 1), mar=c(3.8, 3.6, 2.1, 1.4))

with(family_mean, plot(map_length, log10_size, pch=19, axes=FALSE, cex=0.7,
                       ylab='', xlab='', xlim=c(0, 52), ylim=c(-4.3, 1.1), col='gray48'))

axis(1, seq(0, 50, 10))
axis(2, seq(-4, 1), seq(-4, 1), labels=latex2exp::TeX(sprintf("$10^{%d}$", seq(-4, 1))))

with(bs_p_ci, polygon(c(x, rev(x)), c(lower, rev(upper)), 
                    col=scales::alpha('cornflowerblue', 0.2), lty=0))

with(bs_f_ci, polygon(c(x, rev(x)), c(lower, rev(upper)), 
                    col=scales::alpha('gray42', 0.2), lty=0))

with(phyla_mean, points(map_length, log10_size, pch=19, col='cornflowerblue', cex=0.7))

lines(sf, predict(mf2, newdata=tibble(map_length=sf)), lty='dashed', col='gray32')
lines(sp, predict(mp2, newdata=tibble(map_length=sp)), lty='dashed', col='cornflowerblue')


legend(31, -3.5, inset=0, legend=c('family averages', 'phylum averages'), 
       fill=c('gray42', 'cornflowerblue'),
			 title='',
			 adj=c(0, 0.45),
			 cex=1.2,
			 bty='n', text.col='gray24',
			 border=0)

labels <- with(family_mean, identify(map_length, log10_size, family))

text(12.4, -3.83, "Platyhelminthes", col='cornflowerblue', cex=0.8, adj=0, font=3)
text(4.4, -4.10, "Diplogasteridae", col='gray42', cex=0.8, adj=0, font=3)

title(xlab='map length (Morgans)', ylab='body size (meters)', line=2.5, cex.lab=1.4)

if (pdf) dev.off()
par(opar)

