# map_length_body_size.r -- create the body size / map length figure

library(tidyverse)
library(ggrepel)
library(viridis)
library(grDevices)
library(wesanderson)
library(shape)
library(rstan)
source('../../R/utilities.r')


load('../../data/body_size_map_length_bayesian.Rdata')
opar <- par(no.readonly=TRUE)
arrow_col <- 'gray32'

pdf <- TRUE
# pdf <- FALSE
if (pdf) pdf("map_length_body_size_bayesian.pdf", width=7.3, height=7.3/1.9)

par(mfrow=c(1, 2), mar=c(3.8, 3.6, 2.1, 1.4))

with(family_mean, plot(map_length, log10_size, pch=19, axes=FALSE, cex=0.7,
                       ylab='', xlab='', xlim=c(0, 52), ylim=c(-4.3, 1.1), 
                       type='n', col='gray48'))

lines(x_new, colMeans(extract(bmf3, 'mu_rep')$mu_rep), 
      lty='solid', col=scales::alpha('red', 0.8), lwd=1, lend='butt')

with(family_mean, points(map_length, log10_size, pch=19, cex=0.7,
                        ylab='', xlab='', xlim=c(0, 52), ylim=c(-4.3, 1.1), col='gray48'))


axis(1, seq(0, 50, 10))
axis(2, seq(-4, 1), seq(-4, 1), labels=latex2exp::TeX(sprintf("$10^{%d}$", seq(-4, 1))))

ci_polygon(x_new, bmf1, color='gray42', type='hpdi', 
           transp=0.2, par='mu_rep')
lines(x_new, colMeans(extract(bmf1, 'mu_rep')$mu_rep), 
      lty='dashed', col='gray42')

with(phyla_mean, points(map_length, log10_size, pch=19, col='cornflowerblue', cex=0.7))
ci_polygon(x_new, bmp1, color='cornflowerblue', type='hpdi', 
           transp=0.2, par='mu_rep')
lines(x_new, colMeans(extract(bmp1, 'mu_rep')$mu_rep), 
      lty='dashed', col='cornflowerblue')

# PPD
#ci_polygon(x_new, bmf1, color='gray42', type='quantile', 
#           transp=0.1, par='y_rep')
#ci_polygon(x_new, bmp1, color='cornflowerblue', type='quantile', 
#           transp=0.1, par='y_rep')

legend(27, -3, inset=0, legend=c('family averages', 'phylum averages'), 
       fill=c('gray42', 'cornflowerblue'),
			 title='',
			 adj=c(0, 0.45),
			 cex=0.7,
			 bty='n', text.col='gray24',
			 border=0)

#labels <- with(family_mean, identify(map_length, log10_size, family))

text(12.4, -3.83, "Platyhelminthes", col='cornflowerblue', cex=0.4, adj=0, font=3)
text(4.4, -4.10, "Diplogasteridae", col='gray42', cex=0.4, adj=0, font=3)

title(xlab='map length (Morgans)', ylab='body size (meters)', line=2.5)

mtext("A", 3, at=-3, cex=1.4, xpd=TRUE, font=2)

# second plot
load('../../data/diversity_size_theory.Rdata')

set.seed(1)
#xjitter <- rnorm(length(dtft$size), 0, 0.06)
xjitter <- 0
with(dtft, plot(log10(size) + xjitter, log10(diversity / Ne_N_HK94), 
								type='n', 
                axes=FALSE,
								xlim=c(0.4, -3), 
                ylim=c(-3.2, 0),
								ylab='', xlab=''))

with(dtft, Arrows(log10(size) + xjitter, log10(diversity), 
                  log10(size) + xjitter, log10(diversity / Ne_N_HK94),
									arr.type="triangle", arr.width=0.07, 
                  arr.length=0.05, lcol='gray84',
                  lwd=0.8,
									arr.col='gray84'))

with(dtft, Arrows(log10(size) + xjitter, log10(diversity), 
                  log10(size) + xjitter, log10(diversity / Ne_N_SC98),
									arr.type="triangle", arr.width=0.07, 
                  arr.length=0.05, lcol='gray32',
                  lwd=0.8,
									arr.col=arrow_col))


# abline(a=-2.5, b=0.04, lty=2)

nbins <- 41
ml_cols <- rev(viridis(nbins))
#ml_cols <- rev(wes_palette('Zissou1'))
#ml_cols <- wescols[c(2, 3, 4, 1)]
#mwescolsl_bins <- cut_interval(round(dtft$map_length, 0), 5)
#ml_bins <- cut(round(dtft$map_length, 0), c(0, 10, 20, 30, 41))
ml_bins <- cut_interval(dtft$map_length, n=41)
stopifnot(all(!is.na(ml_bins)))

midpoint <- function(x) {
  out <- lapply(strsplit(x, ','), function(y) {
           as.numeric(gsub('(\\)|\\]|\\[|\\()', '', y))
       })
  rowMeans(do.call(rbind, out))
}

with(dtft, points(log10(size) + xjitter, log10(diversity), 
									pch=19,
									cex=0.8, #col='gray23'))
									col=ml_cols[ml_bins]))

# experimental stuff
# selfers <- c('Caenorhabditis briggsae', 'Caenorhabditis elegans')
# f <- nls(log10(diversity / Ne_N_SC98) ~ a + log10(size), 
#            data=dtft %>% filter(!species %in% selfers),
#            start=list(a=-3))
# x <- seq(0, -3, length.out=100)
# lines(x, -3.2 - x)

# this is ridiculous, but I'm OCD as fuck.... but I want certain points
# to not obscure arrows. 
overdraw <- c(1L, 22L, 26L)
with(dtft[overdraw, ], Arrows(log10(size), log10(diversity), 
                              log10(size), log10(diversity / Ne_N_SC98),
									       arr.type="triangle", arr.width=0.07, arr.length=0.05, 
                         lcol=arrow_col,
                         lwd=0.8,
									       arr.col=arrow_col))
with(dtft[overdraw, ], points(log10(size), log10(diversity), 
								              pch=19,
								              cex=0.8, #col='gray23'))
								              col=ml_cols[ml_bins]))

abline(h=-2, col=alpha('gray42', 0.2), lty=2)
abline(h=-1, col=alpha('gray42', 0.2), lty=2)
abline(h=0, col=alpha('gray42', 0.2), lty=2)
axis(1, seq(0, -3), labels=latex2exp::TeX(sprintf("$10^{%d}$", seq(0, -3))))
axis(2, seq(0, -3), labels=latex2exp::TeX(sprintf("$10^{%d}$", seq(0, -3))))

# a line showing the diversity range
#segments(0.507, log10(min(dtft$diversity)), 0.507, log10(max(dtft$diversity)), 
#         lwd=3, xpd=TRUE, lend='butt',
#         col='gray68')
#
#segments(log10(max(dtft$size)), -3.167, log10(min(dtft$size)), -3.167,
#         lwd=3, xpd=TRUE, lend='butt',
#         col='gray68')

# continuous gradient legend
legend_image <- as.raster(matrix(ml_cols, nrow =1))
lx <- 0.25
ly <- -0.6
width <- 1.75
height <- 0.15
rasterImage(legend_image, lx, ly, lx - width, ly - height)
labs <- round(midpoint(levels(ml_bins)))
# use a subset of colors
i <- labs[seq(1, length(levels(ml_bins)), length.out=4)]
labx <- seq(lx, lx - width, length.out=length(i))
text(labx, ly + 0.1, labs[i], cex=0.6)
text((2*lx - width)/2, ly + 0.25, "map length (Morgans)",
     adj=c(0.5, 0), cex=0.8)

# legend(x=0.5, y=-0.1, inset=0, 
#        legend=labs[i], fill=ml_cols[i],
# 			 title='  map length (Morgans)',
# 			 adj=c(-0.15, 0.5),
#        horiz=TRUE,
# 			 cex=0.7, ncol=1, xpd=TRUE,
# 			 bty='n', text.col='gray24',
# 			 border=0)

#labels <- with(dtft, identify(log10(size), log10(diversity), species))

# labels
text(-2.62, -2.74, "selfing\nCaenor.", cex=0.5, col='gray24', font=3)
text(-2.2,  -1.5, "Anopheles\ngambiae", cex=0.5, col='gray24', font=3)

text(-2.76, -2.14, '{', srt = 90, cex = 0.9, col='gray42')
text(-2.67,  -2.23, "Drosophila", cex=0.5, col='gray24', font=3)
text(-0.15,  -3., "H. sapiens", cex=0.5, col='gray24', font=3)
text(-0.6,  -1.64, "Peromyscus\nmaniculatus", cex=0.5, col='gray24', font=3)
text(-1.26,  -2.8, "Monodelphis\ndomestica", cex=0.5, col='gray24', font=3)
text(-3.1,  -1.64, "Daphnia\nmagna", cex=0.5, col='gray24', font=3, xpd=TRUE)
text(-1.3,  -1.34, "Bombyx\nmandarina", cex=0.5, col='gray24', font=3)
text(-1.8,  -2.3, "Apis\nmellifera", cex=0.5, col='gray24', font=3)

text(-2.5, -3.12, 'large N', adj=0, font=4, cex=0.8, col='gray24', xpd=TRUE)
text(0.4, -3.12, 'small N', adj=0, font=4, cex=0.8, col='gray24')

title(xlab='body size (meters)', ylab='diversity', line=2.5)

mtext("B", 3, at=0.6, cex=1.4, font=2)

if (pdf) dev.off()
par(opar)

