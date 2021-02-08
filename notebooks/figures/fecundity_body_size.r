# fecundity_body_size.r -- create the fecundity / body size figure

library(tidyverse)
library(ggrepel)
library(viridis)
library(wesanderson)
library(shape)
library(rstan)
source('../../R/utilities.r')

arrow_col <- 'gray32'

opar <- par(no.readonly=TRUE)
load('../../data/body_size_fecundity.Rdata')

pdf <- TRUE
#pdf <- FALSE
if (pdf) pdf("fecundity_body_size.pdf", width=7.3, height=7.3/1.9)

par(mfrow=c(1, 2), mar=c(3.8, 3.6, 2.1, 1.4))


with(dr_family_averages %>% filter(fecundity < 8e4),
     plot(log10_size, log10_fecundity, 
          axes=FALSE,
          pch=19,
          xlab='', ylab='',
          xlim=c(-3, 1),
          ylim=c(-3, 5),
          col='gray42'))

with(dr_phylum_averages %>% filter(fecundity < 8e4),
     points(log10_size, log10_fecundity, 
          pch=19,
          col='cornflowerblue'))


# with(dr,
#      points(log10_size, log10_fecundity, 
#           pch=19,
#           cex=0.4,
#           col='gray42'))


axis(1, seq(-3, 1), labels=latex2exp::TeX(sprintf("$10^{%d}$", seq(-3, 1))))
axis(2, seq(-3, 6, 2), labels=latex2exp::TeX(sprintf("$10^{%d}$", seq(-3, 6, 2))))
title(xlab='body size (meters)', ylab='fecundity', line=2.5)

ci_polygon(lnf_xnew, y=extract(lnf, pars='y_rep')$y_rep,
           color='gray42', type='quantile', transp=0.1)
ci_polygon(lnf_xnew, y=extract(lnf, pars='mu_rep')$mu_rep,
           color='gray42', type='hpdi')
lnf_p <- colMeans(as_tibble(extract(lnf, pars='mu_rep')))
lines(lnf_xnew, lnf_p, col=alpha('black', 0.4), lwd=2)

mtext("A", 3, at=-3.22, cex=1.4, xpd=TRUE, font=2)

legend(-3.15, -1.5, inset=0, legend=c('family averages', 'phylum averages'), 
       fill=c('gray42', 'cornflowerblue'),
			 title='',
			 adj=c(0, 0.45),
			 cex=0.7,
			 bty='n', text.col='gray24',
			 border=0)



## Second plot
dr <- read_tsv('../../data/romiguier_et_al_2014_updated.tsv')

set.seed(1)
xjitter <- rnorm(length(dr$size), 0, 0.1)

with(dr, plot(log10(size) + xjitter, log10(4*piS * (2 + fecundity)),
              pch=19, 
              cex=0.7,
              type='n',
              axes=FALSE,
              ylab='', xlab='', 
              xlim=c(0, -3),
              ylim=c(-3, 5),
              col='gray48'))

with(dr, Arrows(log10(size) + xjitter, log10(piS), 
                log10(size) + xjitter, log10(4*piS * (2 + fecundity)),
                arr.type="triangle", arr.width=0.07, arr.length=0.05, 
                lcol=arrow_col,
                lwd=0.8,
                arr.col=arrow_col))
#abline(lm(log10(4*piS * (2 + fecundity)) ~ log10(size), data=dr))
with(dr, points(log10(size) + xjitter, log10(piS),
              pch=21, 
              bg='gray48'))


axis(1, seq(0, -3), line=0.5)
axis(2, seq(-3, 5, 2), labels=latex2exp::TeX(sprintf("$10^{%d}$", seq(-3, 5, 2))), line=0.0)

text(-2.4, -3.3, 'large N', adj=0, font=4, cex=0.8, col='gray24', xpd=TRUE)
text(0, -3.3, 'small N', adj=0, font=4, cex=0.8, col='gray24', xpd=TRUE)


title(xlab='body size (meters)', ylab='diversity', line=2.5)
mtext("B", 3, at=0.18, cex=1.4, font=2)

if (pdf) dev.off()
