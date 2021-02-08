library(tidyverse)
library(wesanderson)
library(RColorBrewer)

range_cat_levels <- c("Island", "Narrow endemic", "Broad endemic", "Cosmopolitan")

d <- read_tsv('../../data/combined_data.tsv') %>%
       filter(is.finite(range), !is.na(range_cat)) %>%
       mutate(log10_range = log10(range)) %>%
       mutate(range_cat = ifelse(range_cat == 'Unknown', NA, range_cat)) %>%
       mutate(range_cat = factor(range_cat, levels=range_cat_levels, ordered=TRUE)) %>%
       arrange(range_cat, range)  %>% 
       # space the groups
       mutate(xpos = 1:nrow(.) + 0.1*as.integer(range_cat)-1)



opar <- par(no.readonly=TRUE)
# pdf('range_descriptions.pdf', width=8, height=5)
pdf('range_categories_long.pdf', height=5, width=8)

par(mar=c(4, 3, 0, 5))
cols <- wes_palette('Darjeeling1')[c(5, 2, 3, 1, 4)]
plot(d$xpos, d$log10_range, pch = 19, col=cols[d$range_cat], 
     xlab = '', ylim=c(2, 8.9),
     ylab = '', axes=FALSE) 
ylabs = seq(2, 8, 2)
axis(2, ylabs, labels=latex2exp::TeX(sprintf("$10^{%d}$", ylabs)),
     cex.axis=0.8, las=1,
     tck=-0.01, hadj=0.5)
text(x=d$xpos, y=2, labels=d$species, 
     col=cols[d$range_cat], xpd=TRUE, 
     las=1,
     srt=90, adj=1, cex=0.4)
legend(70, 4, range_cat_levels, fill=cols, 
       bty='n', box.lwd=0, border=0, cex=1)

title(ylab=latex2exp::TeX('estimated range ($km^2$)'), line=1.7)

# par(mar=c(8, 5, 2, 1))
# cols <- wes_palette('Darjeeling1')[c(5, 2, 3, 1, 4)]
# plot(d$xpos, d$log10_range, pch = 19, col=cols[d$range_cat], 
#      ylab = latex2exp::TeX('estimated range ($km^2$)'),
#      xlab = '', axes=FALSE, ylim=c(2, 9))
# ylabs = seq(2, 8, 2)
# axis(2, ylabs, labels=latex2exp::TeX(sprintf("$10^{%d}$", ylabs)))
# text(x=d$xpos, y=2, labels=d$species, col=cols[d$range_cat], xpd=TRUE, 
#      las=1, srt=90, adj=1, cex=0.4)
# legend(102, 4.2, range_cat_levels, fill=cols, 
#        bty='n', box.lwd=0, border=0, cex=0.7)



dev.off()
par(opar)
