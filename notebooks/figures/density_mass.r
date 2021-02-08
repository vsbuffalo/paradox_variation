library(tidyverse)
library(scales)

opar <- par(no.readonly=TRUE)

pdf("popsize_hists", width=5, height=5)
par(mar=c(5, 5, 2, 1))

load('../../data/post_pred_draws.Rdata')

sps <- c('Homo sapiens', 'Aedes aegypti', 'Drosophila melanogaster',
         'Balaenoptera bonaerensis', 'Canis latrans/lupus', 'Mus musculus')

d <- post_pred_draws %>% 
      mutate(mean_popsize = map_dbl(popsize, mean)) %>% 
      arrange(mean_popsize)
stopifnot(all(sps %in% d$species))

x <- d %>% filter(species %in% sps)

x %>% unnest() %>% ggplot(aes(x=popsize, fill=species)) +geom_histogram()

nsps <- length(unique(x$species))
sps_cols <- hue_pal()(nsps)


hists <- lapply(x$popsize, function(x) hist(x, breaks=seq(0, 24, length.out=100), freq=TRUE))

lapply(seq_along(hists), function(i) {
  add = ifelse(i == 1, FALSE, TRUE)
  plot(hists[[i]], add=add, col=alpha(sps_cols[i], 0.5), 
       border=0, xlim=c(0, 25), axes=FALSE)
})
xseq <- seq(0, 25, 5)
axis(1, xseq, labels=latex2exp::TeX(sprintf("$10^{%d}$", xseq)), cex.axis=0.8, line=0.5)
yseq <- seq(0, 1200, 400)
axis(2, yseq, cex.axis=0.8, line=0.5)


bins <- do.call(rbind, map(hists, 'density'))
barplot(bins, border=0, col=alpha(sps_cols, 1))

dev.off()
