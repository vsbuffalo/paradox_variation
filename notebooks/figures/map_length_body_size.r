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
if (pdf) pdf("map_length_body_size.pdf", width=7.3, height=7.3/1.9)

par(mfrow=c(1, 2), mar=c(3.8, 3.6, 2.1, 1.4))

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


legend(27, -3, inset=0, legend=c('family averages', 'phylum averages'), 
       fill=c('gray42', 'cornflowerblue'),
			 title='',
			 adj=c(0, 0.45),
			 cex=0.7,
			 bty='n', text.col='gray24',
			 border=0)

labels <- with(family_mean, identify(map_length, log10_size, family))

text(12.4, -3.83, "Platyhelminthes", col='cornflowerblue', cex=0.4, adj=0, font=3)
text(4.4, -4.10, "Diplogasteridae", col='gray42', cex=0.4, adj=0, font=3)

title(xlab='map length (Morgans)', ylab='body size (meters)', line=2.5)

mtext("A", 3, at=-3, cex=1.4, xpd=TRUE, font=2)

# second plot
load('../../data/diversity_size_theory.Rdata')
with(dtft, plot(log10(size), log10(diversity / Ne_N_SC98), 
								type='n', 
                axes=FALSE,
                #axes=TRUE,
								xlim=c(0.4, -3), 
                ylim=c(-3.1, -1),
								ylab='', xlab=''))

with(dtft, Arrows(log10(size), log10(diversity), log10(size), log10(diversity / Ne_N_SC98),
									arr.type="triangle", arr.width=0.07, arr.length=0.05, lcol=arrow_col,
                  lwd=0.8,
									arr.col=arrow_col))

# abline(a=-2.5, b=0.04, lty=2)

nbins <- 4
#ml_cols <- (viridis(nbins))
ml_cols <- rev(wes_palette('Zissou1'))
#ml_bins <- cut_interval(round(dtft$map_length, 0), 5)
ml_bins <- cut(round(dtft$map_length, 0), c(0, 10, 20, 30, 41))
stopifnot(all(!is.na(ml_bins)))

with(dtft, points(log10(size), log10(diversity), 
									pch=19,
									cex=0.8, #col='gray23'))
									col=ml_cols[ml_bins]))

# this is ridiculous, but I'm OCD as fuck.... but I want certain points
# to not obscure arrows. 
overdraw <- c(1L, 22L, 26L)
with(dtft[overdraw, ], Arrows(log10(size), log10(diversity), log10(size), log10(diversity / Ne_N_SC98),
									       arr.type="triangle", arr.width=0.07, arr.length=0.05, lcol=arrow_col,
                         lwd=0.8,
									       arr.col=arrow_col))
with(dtft[overdraw, ], points(log10(size), log10(diversity), 
								              pch=19,
								              cex=0.8, #col='gray23'))
								              col=ml_cols[ml_bins]))

axis(1, seq(0, -3), labels=latex2exp::TeX(sprintf("$10^{%d}$", seq(0, -3))))
axis(2, seq(-1, -3), labels=latex2exp::TeX(sprintf("$10^{%d}$", seq(-1, -3))))

# a line showing the diversity range
segments(0.507, log10(min(dtft$diversity)), 0.507, log10(max(dtft$diversity)), 
         lwd=3, xpd=TRUE, lend='butt',
         col='gray68')

segments(log10(max(dtft$size)), -3.167, log10(min(dtft$size)), -3.167,
         lwd=3, xpd=TRUE, lend='butt',
         col='gray68')


legend(x=0.5, y=-0.8, inset=0, legend=levels(ml_bins), fill=ml_cols,
			 title='map length (Morgans)',
			 adj=c(0, 0.45),
			 cex=0.7, ncol=2, xpd=TRUE,
			 bty='n', text.col='gray24',
			 border=0)

labels <- with(dtft, identify(log10(size), log10(diversity), species))

# labels
text(-2.62, -2.7, "selfing\nCaenor.", cex=0.5, col='gray24', font=3)
text(-2.2,  -1.5, "Anopheles\ngambiae", cex=0.5, col='gray24', font=3)

text(-2.76, -2.13, '{', srt = 90, cex = 0.9, col='gray42')
text(-2.67,  -2.21, "Drosophila", cex=0.5, col='gray24', font=3)
text(-0.15,  -3., "H. sapiens", cex=0.5, col='gray24', font=3)
text(-0.6,  -1.7, "Peromyscus\nmaniculatus", cex=0.5, col='gray24', font=3)
text(-1.24,  -2.8, "Monodelphis\ndomestica", cex=0.5, col='gray24', font=3)
text(-3.1,  -1.6, "Daphnia\nmagna", cex=0.5, col='gray24', font=3, xpd=TRUE)
text(-1.3,  -1.4, "Bombyx\nmandarina", cex=0.5, col='gray24', font=3)
text(-1.8,  -2.26, "Apis\nmellifera", cex=0.5, col='gray24', font=3)

text(-2.5, -3.09, 'large N', adj=0, font=4, cex=0.8, col='gray24', xpd=TRUE)
text(0.4, -3.09, 'small N', adj=0, font=4, cex=0.8, col='gray24')

title(xlab='body size (meters)', ylab='diversity', line=2.5)

mtext("B", 3, at=0.6, cex=1.4, font=2)

if (pdf) dev.off()
par(opar)

