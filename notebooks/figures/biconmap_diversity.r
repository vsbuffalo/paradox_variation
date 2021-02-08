library(viridis)
library(phytools)

load('../../data/biconmap_data.Rdata')
source('color_scheme.r')
source('../../R/utilities.r')


# output <- FALSE
output <- TRUE

if (output)
  pdf('biconmap_diversity.pdf', height=8, width=8)

par(oma=c(0, 0, 0.5, 0))
vcols <- viridis(100)
vcols2 <- viridis(100, option='A')
layout(matrix(1:3, 1, 3), widths=c(0.7, 0.16, 0.7))
pscm2 <- setMap(pscm, colors=vcols)
plot(pscm2, fsize=c(0, 0.8), ftype=c("off","reg"), sig=1, lwd=2,
     outline=FALSE, mar=c(1.1,0.1,1.1,0.01), legend=FALSE)
add.color.bar.log(500, vcols, prompt=FALSE, outline=FALSE, 
              subtitle="", digits=0,
              lims=range(dt_dps$log10_popsize), lwd=4,
              title='pop size')
title(main = 'pop size', font.main = 1, cex.main=1.5)
# ylim<-c(1-0.12*(length(tr_mlps$tip.label)-1), length(tr_mlps$tip.label))

plot.new()
plot.window(xlim=c(-0.3,0.3), ylim=c(1, length(tr_dps$tip.label)))
text(rep(0, length(tr_dps$tip.label)), 
     1:length(tr_dps$tip.label), 
     tr_dps$tip.label,
     font=3, cex=0.5, col=phyla_cols[dt_dps$phylum])
dvcm2 <- setMap(dvcm, colors=vcols2)
plot(dvcm2, fsize=c(0, 0.8), ftype=c("off","reg"),
    direction="leftwards", sig=1, lwd=2, legend=FALSE,
    outline=FALSE, mar=c(1.1,0.01,1.1,0.01))
add.color.bar.log(500, vcols2, prompt=FALSE, outline=FALSE, 
              subtitle="", x=250, digits=0,
              lims=range(dt_dps$log10_diversity), lwd=4,
              title='diversity')
title(main = 'diversity', font.main = 1, cex.main=1.5)
nms <- tools:::toTitleCase(unique(dt_dps$phylum)) 
i <- order(nms)
nms <- nms[i]
legend(300, 30, nms,
       fill=phyla_cols[unique(dt_dps$phylum)][i],
       bty='n', border=0, cex=1, ncol=2)

if (output)
  dev.off()



