library(viridis)
library(tidyverse)
library(phytools)

load('../../data/biconmap_data.Rdata')
load('../../data/main_datasets.Rdata')
source('color_scheme.r')
source('../../R/utilities.r')

# output <- FALSE
output <- TRUE
if (output) 
  pdf('popsize_tree.pdf', height=5, width=4)

par(oma=c(0.2, 0.2, 0.2, 0.2), mar=c(1, 1, 1, 1))
vcols <- viridis(100)
vcols2 <- viridis(100, option='A')
layout(matrix(c(1, 2, 3, 1, 2, 3), byrow=TRUE, ncol=3, nrow=2),  
       width=c(1, 0.05, 1))
pscm2 <- setMap(pscm, colors=vcols)
plot(pscm2, fsize=c(0, 0.8), ftype=c("off","reg"), sig=1, lwd=1.4,
     outline=FALSE, mar=c(1.1,0.1,1.1,0.01), legend=FALSE)
add.color.bar.log(500, vcols, prompt=FALSE, outline=FALSE, 
              subtitle="",
              x=10, digits=0,
              lims=range(dt_dps$log10_popsize), lwd=4,
              title='pop size')
title(main = 'pop size', font.main = 1, cex.main=1.5)

plot.new()
width <- 0.5
plot.window(xlim=c(0, width), ylim=c(1, length(tr_dps$tip.label)))

# this is dumb... it shouldn't be this hard, but it is.
tmp <- tibble(species = tr_dps$tip.label) %>%
  left_join(dt_dps %>% select(species, phylum))  

# get rectangle coords
tmp2 <- tmp %>%
  group_by(phylum) %>% summarize(n=n()) %>% mutate(left=0) 

# order by tip phyla
# get in the same order as data (grrr...)
order_phy <- unique(tmp$phylum)
tmp2 <- tmp2[match(order_phy, tmp2$phylum), ] 

tmp2 <- tmp2 %>% mutate(cumsum = cumsum(n), left = lag(cumsum, default=0)) 
rect_z <- cbind(left=tmp2$left, right=tmp2$cumsum)
rownames(rect_z) <- tmp2$phylum


offset <- 0.5
for (i in 1:nrow(rect_z)) {
  rect(0, offset + rect_z[i, 1], width, offset + rect_z[i, 2],
       lwd=0, border=NA,
       col=alpha(phyla_cols[rownames(rect_z)[i]], 0.9))
  print(phyla_cols[rownames(rect_z)[i]])
}

dvcm2 <- setMap(dvcm, colors=vcols2)
plot(dvcm2, fsize=c(0, 0.8), ftype=c("off","reg"),
    direction="leftwards", sig=1, lwd=1.4, legend=FALSE,
    outline=FALSE, mar=c(1.1,0.01,1.1,0.01))
add.color.bar.log(500, vcols2, prompt=FALSE, outline=FALSE, 
              subtitle="", x=250,
              lims=range(dt_dps$log10_diversity), lwd=4,
              title='diversity')
title(main = 'diversity', font.main = 1, cex.main=1.5)
nms <- tools:::toTitleCase(unique(dt_dps$phylum)) 
i <- order(nms)
nms <- nms[i]
legend(200, 33, nms,
       fill=phyla_cols[tmp2$phylum][i],
       bty='n', border=0, cex=0.8, ncol=2)


if (output) 
  dev.off()

