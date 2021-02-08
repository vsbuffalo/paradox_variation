library(viridis)
library(phytools)

load('../../data/biconmap_maplen_data.Rdata')

source('color_scheme.r')
source('../../R/utilities.r')

pdf('biconmap.pdf', height=8, width=8)
par(oma=c(0, 0, 0.5, 0))
vcols <- viridis(100)
vcols2 <- viridis(100, option='A')
layout(matrix(1:3, 1, 3), widths=c(0.7, 0.14, 0.7))
pscm <- psmlcm
pscm2 <- setMap(pscm, colors=vcols)
plot(pscm2, fsize=c(0, 0.8), ftype=c("off","reg"), sig=1, lwd=2,
     outline=FALSE, mar=c(1.1,0.1,1.1,0.01), legend=FALSE)
add.color.bar.log(500, vcols, prompt=FALSE, outline=FALSE, 
              subtitle="", digits=0,
              lims=range(dt_ml$log10_popsize), lwd=4,
              title='log10(pop size)')
title(main = 'log10(pop size)')
# ylim<-c(1-0.12*(length(tr_mlps$tip.label)-1), length(tr_mlps$tip.label))

plot.new()
plot.window(xlim=c(-0.3,0.3), ylim=c(1, length(tr_ml$tip.label)))
text(rep(0, length(tr_ml$tip.label)), 1:length(tr_ml$tip.label), 
     tr_ml$tip.label,
     font=3, cex=0.5, col=phyla_cols[factor(dt_ml$phylum)])
mlcm2 <- setMap(mlcm, colors=rev(vcols2))
plot(mlcm2, fsize=c(0, 0.8), ftype=c("off","reg"),
    direction="leftwards", sig=1, lwd=2, legend=FALSE,
    outline=FALSE, mar=c(1.1,0.01,1.1,0.01))
add.color.bar.log(500, rev(vcols2), prompt=FALSE, outline=FALSE, 
              subtitle="", x=250, digits=0,
              lims=range(dt_ml$log10_map_length), lwd=4,
              title='log10(map length)')
title(main = 'log10(map length)')
legend(350, 15, tools:::toTitleCase(unique(dt_ml$phylum)), 
       fill=phyla_cols,
       bty='n', border=0, cex=1, ncol=2)

dev.off()



