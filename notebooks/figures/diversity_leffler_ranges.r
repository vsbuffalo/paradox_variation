library(tidyverse)
library(wesanderson)
library(RColorBrewer)

range_cat_levels <- c("Island", "Narrow endemic", "Broad endemic", "Cosmopolitan")

d <- read_tsv('../../data/combined_data.tsv') %>%
       filter(is.finite(log10_diversity), is.finite(range), !is.na(range_cat)) %>%
       mutate(log10_range = log10(range)) %>%
       mutate(range_cat = ifelse(range_cat == 'Unknown', NA, range_cat)) %>%
       mutate(range_cat = factor(range_cat, levels=range_cat_levels, ordered=TRUE)) %>%
       arrange(range_cat, log10_diversity) %>%
       # space the groups
       mutate(xpos = seq(1, nrow(.)) + 2.2*as.integer(range_cat)-1) 

means <- d %>% group_by(range_cat) %>% 
  summarize(xmin=min(xpos), 
            xmax=max(xpos),
            log10_diversity = mean(log10_diversity, na.rm=TRUE))

# ggplot(d, aes(xpos, log10_diversity, color=range_cat, size=10^log10_density)) + geom_point()  + facet_wrap(~range_cat)  + geom_text_repel(aes(xpos, log10_diversity, label=species))

opar <- par(no.readonly=TRUE)
# pdf('range_descriptions.pdf', width=8, height=5)
pdf('diversity_leffler_ranges.pdf', width=5, height=8)

par(mar=c(4, 1, 1, 5))
cols <- wes_palette('Darjeeling1')[c(5, 2, 3, 1, 4)]

scaler <- function(x, range=c(1, 3)) {
  rng <- range(x, na.rm=TRUE)
  slope = abs(diff(range)) / (rng[2] - rng[1])
  min(range) + slope*x
}

plot(d$log10_diversity, d$xpos, pch = 21, col='white',
     bg=cols[d$range_cat], 
     lwd=0.7,
     axes=FALSE,
     cex=scaler(d$log10_density),
     xlab = '',
     ylab = '')

with(means, segments(y0=xmin, y1=xmax, 
                     x0=log10_diversity, x1=log10_diversity,
                     col=cols[range_cat], lwd=1, lty=1))

ylabs = seq(-4, -1)
axis(1, ylabs, labels=latex2exp::TeX(sprintf("$10^{%d}$", ylabs)),
     cex.axis=0.8,
     line=0.5,
     tck=-0.01, padj=-1)
text(y=d$xpos, x=-1, labels=d$species, 
     col=cols[d$range_cat], xpd=TRUE, 
     las=1, srt=0, adj=0, cex=0.3)
legend(-4.2, 55, range_cat_levels, fill=cols, 
       bty='n', box.lwd=0, border=0, cex=1)

sizes <- c(-1, 0, 1, 3, 5, 7)
pex <- round(scaler(sizes), 1)
legend(-4.1, 85, latex2exp::TeX(sprintf("$10^{%0.f}$", sizes)), 
       title.adj=1,
       pt.cex=pex, pch=19, y.intersp=1.2,
       bty='n', box.lwd=0, border=0, cex=1)

text(-3.5, 85, latex2exp::TeX("individuals per km$^2$"), cex=1)

title(xlab=latex2exp::TeX('pairwise diversity (differences per bp)'), line=2.1)

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
