library(viridis)
library(tidyverse)
library(phytools)
library(scales)

# note: why not two-tail? Here I don't think it's appropriate since
# I use the abs. value of the contrasts

load('../../data/nodeheights.Rdata')
# for EDA at end (to dig around and find what the outlier nodes are)
load('../../data/main_datasets.Rdata')
source('color_scheme.r')
source('../../R/utilities.r')
source('../../R/node_height.r')

# output <- FALSE
output <- TRUE
if (output) 
  pdf('node_heights.pdf', height=3, width=4)



layout(matrix(c(1, 2, 3, 
                1, 4, 5, 
                1, 6, 6), byrow=TRUE, ncol=3, nrow=3), 
       width=c(0.25, 1, 1), height=c(1, 1, 0.25))

par(mar=c(1.3, 1.3, 1.2, 1), oma=c(0.1, 0.1, 0.1, 0.1))
plot.new()
title(ylab='absolute value of contrasts', line=-1)

plot_nh(dvnh, 'diversity', ylim=c(0, 0.006), 
        xlab='', ylab='')
plot_nh(psnh, 'pop size', ylim=c(0, 1.1), 
        ylab='', 
        xlab='')
plot_nh(bmnh, 'body mass', xlab='', ylab='', 
        ylim=c(0, 0.4))
plot_nh(rnnh, 'range', ylim=c(0, 1), 
        ylab='', 
        xlab='')

# what's excluded?
dvnh$d  %>% filter(pic > 0.01) %>% pull(tips)
dvnh$d  %>% filter(pic > 0.002) %>% pull(tips)
psnh$d  %>% filter(pic > 2) %>% pull(tips)
rnnh$d  %>% filter(pic > 1) %>% pull(tips)
tips <- bmnh$d  %>% filter(pic > 0.5) %>% pull(tips)

# these all look reasonable -- just lots of variability
dt_dps %>% filter(species %in% tips[[1]]) %>%
  select(species, size, log10_body_mass)


# extra eda
dvnh$d  %>% filter(pic > 0.002) %>% pull(tips)

# now, looking right to left, what are the outliers in the 
# diversity node height test?

# recent outliers
dvnh$d  %>% filter(pic > 0.001, abs(bt - 100) < 50) %>% pull(tips)

# outliers, 200 Mya
tips <- dvnh$d  %>% filter(pic > 0.0015, abs(bt - 200) < 100) %>% pull(tips)

dt_dps %>% filter(species %in% tips[[2]]) %>%
  select(species, size, log10_popsize, diversity, 
         log10_body_mass)

dvnh$d  %>% filter(pic > 0.0015, abs(bt - 200) < 100) %>% pull(bt)
# list(c("Ciona intestinalis", "Ciona savignyi"), 
# c("Bostrycapulus aculeatus", "Crepidula plana", "Crepidula fornicata"))

# far outliers, ~ 570 Mya
dvnh$d  %>% filter(pic > 0.001, abs(bt - 550) < 100) 
tips <- dvnh$d  %>% filter(pic > 0.001, abs(bt - 550) < 100) %>% pull(tips)
unique(dt_dps$class[dt_dps$species %in% tips[[1]]])
# this is the Ascidiacea - non-Ascidiacea split.

unique(dt_dps$class[dt_dps$species %in% tips[[2]]])
# this is cephlapoda - non-cephlaoda split.

bts <- dvnh$d  %>% filter(pic > 0.001, abs(bt - 550) < 100) %>% pull(bt)
dt_dps$class[dvnh$d$bt  == bts]

# body mass outliers
tips <- bmnh$d  %>% filter(pic > 0.6) %>% pull(tips)
dt_dps$family[dt_dps$species %in% tips[[1]]]

tips <- bmnh$d  %>% filter(pic > 0.4, pic< 1) %>% pull(tips)
dt_dps$genus[dt_dps$species %in% tips[[1]]]


# plot_nh(dvnh, 'diversity', ylim=c(0, 0.008))

plot.new()
title(xlab='branching time (Mya)', line=-0.8, xpd=TRUE)

if (output) 
  dev.off()
