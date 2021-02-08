library(scales)
library(tidyverse)
library(shape)
library(ggrepel)


opar <- par(no.readonly=TRUE)

source('color_scheme.r')
source('../../R/utilities.r')

dests <- read_tsv('../../data/sweep_param_ests.tsv')

dshap_raw <- read_tsv('../../data/shapiro_et_al_2007.tsv', col_names=FALSE)

# unsure what all the columns are, but this is what matters
dshap <- tibble(rec = dshap_raw$X3, pi = dshap_raw$X4)


hh_nls <- nls(pi ~ neutral.pi*rec/(rec + alpha), dshap, 
              start=list(neutral.pi=2, alpha=1e-9))

pdf <- TRUE
# pdf <- FALSE
if (pdf) pdf("linked_sel_params.pdf", width=7.3, height=7.3/2.3)

par(mfrow=c(1, 2), mar=c(1, 1, 1, 1), oma=c(1, 1, 1, 1))

plot(log10(dests$nu), rep(1, nrow(dests)), xlim=c(-12, -8), 
    ylab='', xlab='', pch=19, axes=FALSE, 
    col=ifelse(dests$study == '***', 'red', 'black'))

xlabs <- seq(-12, -8)
axis(1, xlabs, line=-5.95,
     cex.axis=0.7,
     tck=-0.02,
     padj=-1.2,
     labels=latex2exp::TeX(sprintf("$10^{%d}$", xlabs)))
mtext(latex2exp::TeX("sweeps per basepair, per generation ($\\nu_{BP})"),
      1, line=-3.4)
mtext("A", 3, at=-12, line=-0, cex=1.4, font=2)

fs <- 0.5
lwd <- 0.8
lcol <- 'gray42'
text(-11.5, 1.05, "Macpherson et al. (2007)", cex=fs)
segments(-11.5, 1.03, -11.45, 1, lwd=lwd, col=lcol)
# text(-10.8, 1.1, "***", cex=fs)
# segments(-10.8, 1.08, -10.61, 1, lwd=lwd, col=lcol)
text(-9.8, 1.05, "Jensen et al. (2008)", cex=fs)
segments(-10.4, 1.08, -10.4, 1, lwd=lwd, col=lcol)
text(-10, 1.1, "Li and Stephan (2006)", cex=fs)
segments(-10.1, 1.03, -10.13, 1, lwd=lwd, col=lcol)
text(-8.5, 1.1, "Andolfatto (2007)", cex=fs)
segments(-8.5, 1.08, -9.15, 1, lwd=lwd, col=lcol)

# overlay
points(log10(dests$nu), rep(1, nrow(dests)), 
    pch=19, col=ifelse(dests$study == '***', 'red', 'black'))


# dests %>%
#   mutate(color=ifelse(study == '***', 'red', 'black')) %>%  
#   ggplot() + 
#   geom_hline(yintercept=0, color='gray32') + 
#   geom_point(aes(nu, 0, color=color), size=4) + 
#   geom_text_repel(aes(nu, 0, label=study), point.padding=1, box.padding=3, seed=1) +
#   scale_x_log10() +   scale_colour_identity()+
#   xlab("beneficial substitutions per bp (ν_bp)") + ylab('')



par(mar=c(3, 3, 1, 0))

plot(pi ~ rec, dshap, pch=16, xlab='', 
     cex=0.7,
     col=scales::alpha('gray21', 0.4),
     axes=FALSE,
     ylab='')
title(ylab=latex2exp::TeX('$\\pi$'), xlab='recombination (cM/Mb)', line=2)
rec <- seq(0, 2.5, length.out=100)
y <- predict(hh_nls, data.frame(rec=rec))
xlab <- seq(0, 2.5, 0.5)
axis(1, xlab, xlab, 
     padj = -1.3,
     line=0.2,
     cex.axis=0.8,
     tck=-0.02)

ylab <- seq(0, 0.06, 0.02)
axis(2, ylab, 
     las=1, 
     cex.axis=0.8,
     line=0.1,
     tck=-0.02, hadj=0.65,
     labels=ylab)

pi0 <- coef(hh_nls)['neutral.pi']
alpha <- coef(hh_nls)['alpha']
lines(rec, pi0/(1 + alpha/rec), col='blue', lwd=3, lend='butt')

rs <- 0.92 # from Elyashiv et al. supp
lines(rec, pi0/(1 + rs/rec), col='green', lwd=3, lty=2, lend='butt')

legend(0, 0.06, c('non-linear least squares', 
                  # latex2exp::TeX("My $J_{2,2}$"),
                  latex2exp::TeX('Elyashiv et al. $r_S$')),
        fill=c('blue', 'green'),
        bty='n', border=0, cex=0.6, ncol=1)
mtext("B", 3, at=0, line=-0, cex=1.4, font=2)

dev.off()
