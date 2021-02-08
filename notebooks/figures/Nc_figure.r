# Nc_figure.r

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

# for stan fit
load('../../data/diversity_popsize_chains.Rdata')
load('../../data/main_datasets.Rdata')
phi <- mean(rstan:::extract(dps_fit, pars='phi')$phi)
beta <- mean(rstan:::extract(dps_fit, pars='beta')$beta)

# alias this -- less typing
d <- da_dps %>% filter(!is.na(log10_popsize)) %>% 
      group_by(phylum) %>%  
      mutate(n = n(), phyla_mean = mean(log10_popsize)) %>%
      filter(n > 5) %>%
      arrange(phyla_mean, log10_popsize)

d %>% ggplot() + geom_histogram(aes(log10_popsize, fill=phylum), binwidth=0.7, alpha=0.5)


d %>% ggplot() + geom_density(aes(log10_popsize, fill=phylum), binwidth=0.7, alpha=0.5)

# what to include
# x <- d$pred_log10_N  # the simplistic unregularized version
lx <- 4
ux <- 18


output <- TRUE
# output <- FALSE

if (output)
  pdf("Nc_figure.pdf", height=7, width=7/1.5)


par(mar=c(1, 1, 1, 1), oma=c(0, 0, 0, 0))



nf <- layout(matrix(c(1, 2), ncol=1), heights=c(2, 8))
x <- d$log10_popsize
sps <- d$species
phy <- d$phylum
i <- seq_along(x)

## first plot
bins = seq(lx, ux, length.out=14)
hd <- d %>% nest() %>% 
  mutate(hist = map(data, ~ hist(.$log10_popsize, bins, plot=FALSE)))

al <- 1

par(mar=c(1, 1, 1, 1))

h <- hist(d$log10_popsize, 10, plot=FALSE)
xd <- do.call(rbind, map(hd$hist, 'counts'))
idx <- c(3, 2, 1, 4)
xd <- xd[idx, ]
cols <- phyla_cols[hd$phylum[idx]]

barplot(xd, axes=FALSE, horiz=FALSE,
        space=0, border='white', col=cols)

plot(x, i, col=phyla_cols[d$phylum], pch=19, 
     ylab='',
     xlab='',
     axes=FALSE,
     xlim=c(lx, ux))
title(xlab='approximate population size', line=-0.4, cex.lab=1.2)

xseq <- seq(lx, ux, 2)

labs <- latex2exp::TeX(sprintf("$10^{%d}$", xseq)) 
# this is some dumb shit to align text...
labs[1:3] <- latex2exp::TeX(sprintf("  $10^{%d}$", xseq[1:3])) 
axis(3, xseq, labels=labs, 
     line=1.4,  padj=3,
     las=1, tck=0.01)
axis(3, xseq, labels=rep("", 8),
     tck=-0.01,
     line=1.4)


# labels
keep_sps <- c("Homo sapiens", 'Drosophila melanogaster',
              'Gorilla beringei/gorilla', 'Eschrichtius robustus',
              'Pan paniscus', 'Lynx lynx',
              'Gulo gulo', 'Ciona savignyi', 'Drosophila sechellia', 
              'Culex pipiens', 'Aedes aegypti', 
              'Halictus scabiosae', 'Pheidole pallidula',
              'Drosophila arizonae', 'Mytilus edulis', 
              'Heliconius melpomene',
              'Crassostrea gigas', 'Mytilus californianus', 'Strongylocentrotus purpuratus',
              'Echinocardium cordatum', 'Apis mellifera', 'Cystodytes dellechiajei')
labs <- c("Homo sapiens", 'Drosophila melanogaster',
          'Gorilla gorilla', 'Eschrichtius robustus',
          'Pan paniscus', 'Lynx lynx',
          'Gulo gulo', 'Ciona savignyi', 'Drosophila sechellia', 
          'Culex pipiens', 'Aedes aegypti', 
           'Halictus scabiosae', 'Pheidole pallidula',
           'Drosophila arizonae', 'Mytilus edulis', 
            'Heliconius melpomene',
          'Crassostrea gigas', 'Mytilus californianus', 'Strongylocentrotus purpuratus',
          'Echinocardium cordatum', 'Apis mellifera', 'Cystodytes dellechiajei')



dput(rep(c(2, 4), length(keep_sps)/2))
sides <- c(4, 4, 2, 4, 2, 4, 2, 1, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 3)
idx <- match(keep_sps, sps) 
# keep_sps[is.na(idx)]
xs <- i[idx]
ys <- d$log10_popsize[idx]
or <- order(xs)
xs <- xs[or]
ys <- ys[or]
labs <- labs[or]
text(ys, xs, labs, srt=0, pos=sides, cex=0.5, xpd=TRUE)

phyla_cols_s <- phyla_cols[unique(d$phylum)]
legend(14, 32, names(phyla_cols_s), fill=phyla_cols_s,
       bty='n', border=0, cex=0.7, ncol=1)



if (output) dev.off()
par(opar)
