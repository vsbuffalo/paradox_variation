library(viridis)
library(grDevices)
library(wesanderson)
library(shape)
library(rstan)

opar <- par(no.readonly=TRUE)

# source('../../R/utilities.r')
source('color_scheme.r')

load('../../data/corbett_detig_data.Rdata')
load('../../data/linked_sel_div_changes.Rdata')

pdf <- TRUE
# pdf <- FALSE
if (pdf) pdf("corbett_detig.pdf", width=7.3, height=7.3/1.9)

par(mfrow=c(1, 2), mar=c(3, 3, 2, 1.1))
# nf <- layout(matrix(c(1, 2), ncol=2))
# layout.show(nf)

plot(log10_diversity ~ log10_popsize, dca, pch=19, axes=FALSE, 
     type='n',
     ylim=c(-3, -1), xlim=c(5.5, 16))
cd_div_fit <- lm(log10(diversity) ~ log10_popsize, dca_noself)
abline(cd_div_fit,  col='gray42')
points(log10_diversity ~ log10_popsize, dca, pch=19, 
     col=phyla_cols[dca$phylum])
phyla_cols_s <- phyla_cols[unique(dca_noself$phylum)]
legend(8, -1, names(phyla_cols_s), 
       fill=phyla_cols_s, bty='n', border=0, cex=0.6, ncol=2)

xseq <- seq(6, 18, 2)
axis(1, xseq, 
     labels=latex2exp::TeX(sprintf("$10^{%d}$", xseq)), 
     padj = -0.9,
     tck=-0.02,
     cex.axis=0.8, line=0.3)
yseq <- seq(-3, -1)
axis(2, yseq, labels=yseq, cex.axis=0.8, line=0, 
     tck=-0.02, hadj=0.2,
     las=1)
title(ylab="diversity",
      xlab="approximate population size", line=1.8)

plot(dml_full$log10_popsize, log10(dml_full$Ne_N_RHH_BGS), 
     axes=FALSE, pch=19,
     col=wes_palette("Darjeeling1")[3],
     xlim=c(6, 18))

points(dca_noself$log10_popsize, log10(1-dca_noself$impact_of_sel), pch=19, 
       col=wes_palette("Darjeeling1")[1])
points(dml$log10_popsize, log10(dml$Ne_N_RHH_BGS), pch=19,
       col=wes_palette("Darjeeling1")[2])

axis(1, xseq, 
     labels=latex2exp::TeX(sprintf("$10^{%d}$", xseq)), 
     padj = -0.9,
     tck=-0.02,
     cex.axis=0.8, line=0.3)
yseq <- seq(-12, 0, 4)
axis(2, yseq, 
     labels=latex2exp::TeX(sprintf("$10^{%d}$", yseq)), 
     cex.axis=0.8, line=0, 
     tck=-0.02, hadj=0.55,
     las=1)
# ios_fit <- lm(impact_of_sel ~ log10_popsize, dca_noself)
# abline(ios_fit, col='gray42')
# points(log10(impact_of_sel) ~ log10_popsize, dca_noself, pch=19, 
#      col=phyla_cols[dca$phylum])

cols = wes_palette('Darjeeling1')[1:3]
legend(6, -10, 
       c( latex2exp::TeX("Corbett-Detig et al."),
         latex2exp::TeX("RHH + BGS, $N = \\pi / 4\\mu$"),
         latex2exp::TeX("RHH + BGS, $N = N_c$")),
       fill=cols, bty='n', border=0, cex=0.6, ncol=1)


title(ylab=latex2exp::TeX("reduction in diversity, $R = N_e / N$"),
      xlab="approximate population size", line=1.8)


#plot(map_length ~ log10_popsize, dat_ws, type='n', 
#     axes=FALSE,
#     ylim=c(0, 60),
#     xlim=c(6, 16.5), 
#     ylab="", 
#     xlab="")
#ci_polygon(x, fit=ln_fit, col='cornflowerblue')
#ci_polygon(x, fit=ln_fit, par='mu_rep', color='gray32')
#polygon(x=c(5, 8, 8, 5), y=c(60, 60, 63, 63),
#        col='white', border=NA)
#points(map_length ~ log10_popsize, dat_ws,
#       col='gray12', cex=0.8, lwd=1,
#       bg='gray23', pch = c(21, 24)[factor(social)])
#lines(log10_popsize_seq, pred)

#alpha <- 1.5
## L = -α / (log(10) * (φ + (β1-1) * log10(N)))
## singularity when 
## 
## lines(x, -28  / (log(10)*(4 - x)), col='red')
#b0 <- coef(div_fit)[1]
#b1 <- coef(div_fit)[2]
#mu <- 10^(-9)
## sapply(c(10^(-8:-10)), function(mu) {
#  # phi <- b0 - log10(mu)
#  # sing <- - phi / (b1 - 1)
#  # x <- seq(sing, 17, length.out=10000)
#  # lines(x, -alpha / (log(10) * (phi + (b1-1)*x)), lwd=1.4,
#        # lty = 2, col='red')
## })
#phi <- b0 - log10(mu)
#sing <- - phi / (b1 - 1)
#x <- seq(sing, 17, length.out=10000)
#lines(x, -alpha / (log(10) * (phi + (b1-1)*x)), lwd=1.4,
#      lty = 2, col='red')
#title(ylab="recombination map length (Morgans)",
#      xlab="approximate population size", line=1.8)

#xseq <- seq(6, 16, 2)
#axis(1, xseq, 
#     labels=latex2exp::TeX(sprintf("$10^{%d}$", xseq)), 
#     padj = -0.9,
#     tck=-0.02,
#     cex.axis=0.8, line=0.3)
#yseq <- seq(0, 60, 10)
#axis(2, yseq, labels=yseq, cex.axis=0.8, line=0, 
#     tck=-0.02, hadj=0.2,
#     las=1)
#mtext("A", 3, at=5.5, line=-0, cex=1.4, font=2)

## text(13, 30, "social insects", adj=c(0, NA), cex=0.6)
## legend(13, 55, "social insects", pch=24, 
## bty='n', border=0, cex=0.8, ncol=2)
#par(mar=c(3, 2.8, 2, 0.1))


### subfigure B
#set.seed(1)
#xjitter <- 0
#arrow_col <- 'gray32'
#with(dat_ws, plot(log10_popsize + xjitter,
#                log10_SC98_diversity, 
#								type='n', 
#                axes=FALSE,
#								xlim=c(5, 17.5), 
#                ylim=c(-3.2, -0.9),
#								ylab='', 
#                xlab=''))

#title(xlab="approximate population size", line=1.8)
#title(ylab = "diversity", line=1.9)


#with(dat_ws, Arrows(log10_popsize + xjitter, 
#                    log10_diversity, 
#                  log10_popsize + xjitter, 
#                  log10_HK94_diversity,
#									arr.type="triangle", arr.width=0.07, 
#                  arr.length=0.05, lcol='gray84',
#                  lwd=0.8,
#									arr.col='gray84'))

#with(dat_ws, Arrows(log10_popsize + xjitter, log10_diversity, 
#                  log10_popsize + xjitter, log10_SC98_diversity,
#									arr.type="triangle", arr.width=0.07, 
#                  arr.length=0.05, lcol='gray32',
#                  lwd=0.8,
#									arr.col=arrow_col))

#nbins <- 40
#ml_bins <- cut_interval(dat_ws$map_length, n=nbins)
#stopifnot(all(!is.na(ml_bins)))


#ml_cols <- rev(viridis(nbins, option='C'))
#with(dat_ws, points(log10_popsize + xjitter, 
#                    log10(diversity), 
#   									pch=19,
#   									cex=0.8, #col='gray23'))
#   									col=ml_cols[ml_bins]))

## this is ridiculous, but I'm OCD as fuck.... but I want certain points
## to not obscure arrows. 
## overdraw <- c(1L, 22L, 26L)
## with(dtft[overdraw, ], Arrows(log10(size), log10(diversity), 
##                               log10(size), log10(diversity / Ne_N_SC98),
## 									       arr.type="triangle", arr.width=0.07, arr.length=0.05, 
##                          lcol=arrow_col,
##                          lwd=0.8,
## 									       arr.col=arrow_col))
## with(dtft[overdraw, ], points(log10(size), log10(diversity), 
## 								              pch=19,
## 								              cex=0.8, #col='gray23'))
## 								              col=ml_cols[ml_bins]))

## abline(h=-2, col=alpha('gray42', 0.2), lty=2)
## abline(h=-1, col=alpha('gray42', 0.2), lty=2)
## abline(h=0, col=alpha('gray42', 0.2), lty=2)
#axis(1, seq(6, 16, 2), line=0.3,
#     padj = -0.9,
#     cex.axis=0.8,
#     tck=-0.02,
#     labels=latex2exp::TeX(sprintf("$10^{%d}$", seq(6, 16, 2))))
#axis(2, seq(-1, -3), 
#     las=1, 
#     cex.axis=0.8,
#     tck=-0.02, hadj=0.65,
#     labels=latex2exp::TeX(sprintf("$10^{%d}$", seq(-1, -3))),
#     line=-0.5)

#legend_image <- as.raster(matrix(ml_cols, nrow =1))
#lx <- 12 
#ly <- -3
#width <- 5
#height <- 0.1
#rasterImage(legend_image, lx, ly, lx + width, ly - height)
#labs <- round(midpoint(levels(ml_bins)))

## use a subset of colors
#i <- seq(1, length(labs), length.out=4)
#labx <- seq(lx, lx + width, length.out=length(i))
## with this number of bins, the labels are 2, 21, 40, and 60
## I round these out to cleaner numbers here:
## text(labx, ly - 0.23, labs[i], cex=0.6)
#text(labx, ly - 0.15, c(2, 20, 40, 60), cex=0.6)
#text((2*lx + width)/2, ly + 0.05, "map length (Morgans)",
#     adj=c(0.5, 0), cex=0.7)
# # ggplot(dat_ws, aes(log10_popsize, log10_diversity)) + geom_point() + geom_text_repel(aes(log10_popsize, log10_diversity, label=species)) 

#lab_cex <- 0.5
#text(16.25, -2.62, "C. elegans", cex=lab_cex, col='gray24', font=3)

#text(14.3,  -1.64, "Anopheles\ngambiae", cex=lab_cex, col='gray24', font=3)
#text(12.99,  -1.46, "Culex pipiens", cex=lab_cex, col='gray24', font=3)
#text(9,  -1.3, "Magallana gigas", cex=lab_cex, col='gray24', font=3)
#text(7.75,  -1.61, "Ciona intestinalis", cex=lab_cex, col='gray24', font=3)

## text(-2.76, -2.14, '{', srt = 90, cex = 0.9, col='gray42')
#text(16.5,  -1.93, "D. melanogaster", cex=lab_cex, col='gray24', font=3)
#text(8,  -3.05, "H. sapiens", cex=lab_cex, col='gray24', font=3)
#text(6,  -3, "Bos taurus", cex=lab_cex, col='gray24', font=3)
## text(-0.6,  -1.64, "Peromyscus\nmaniculatus", cex=lab_cex, col='gray24', font=3)

#text(9.9,  -2.9, "Monodelphis\ndomestica", cex=lab_cex, col='gray24', font=3)
#text(11,  -1.4, "Bombyx\nmandarina", cex=lab_cex, col='gray24', font=3)
#text(13,  -2.25, "Apis\nmellifera", cex=lab_cex, col='gray24', font=3)
#mtext("B", 3, at=4.8, line=-0., cex=1.4, font=2)

if (pdf) dev.off()
par(opar)
 
