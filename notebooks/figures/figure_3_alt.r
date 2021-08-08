library(viridis)
library(ggrepel)
library(grDevices)
library(wesanderson)
library(scales)
library(tidyverse)
library(shape)

opar <- par(no.readonly=TRUE)

source('color_scheme.r')
source('../../R/utilities.r')

load('../../data/linked_sel_div_changes.Rdata')
load('../../data/lognorm_popsize_chains.Rdata')
load('../../data/diversity_popsize_chains.Rdata')
load('../../data/main_datasets.Rdata')


midpoint <- function(x) {
  out <- lapply(strsplit(x, ','), function(y) {
           as.numeric(gsub('(\\)|\\]|\\[|\\()', '', y))
       })
  rowMeans(do.call(rbind, out))
}

pdf <- TRUE
# pdf <- FALSE
if (pdf) pdf("figure_3_alt.pdf", width=7.3, height=7.3/2.3)

par(mfrow=c(1, 2), mar=c(3, 3, 2, 1), oma=c(0, 0, 0, 1))

## subfigure A


# I = 1 - pibar/p0
# pibar / p0 = I - 1
dml %>% filter(species == 'Drosophila melanogaster') %>%
  select(impact_of_sel, Ne_N_CD15) %>% 
  mutate(red = 1-impact_of_sel)

# # if we wantd to use the Elyashiv et al data
# dmel <- 'Drosophila melanogaster'
# dml[dml$species == dmel, ]$Ne_N_CD15 <- 0.77

# dml[dml$species == dmel, ]$log10_CD15_diversity <- log10(dml$diversity[dml$species == dmel]  / dml$Ne_N_CD15[dml$species == dmel])

four_alleles <- function(N, mu) {
  theta <- 4*N*mu
  theta / (1 + (4*theta)/3)
}

dml <- dml_full %>% 
  mutate(Ne_RHH_BGS_strongsel = 10^log10_popsize * Ne_N_RHH_BGS_strongsel) %>%
  mutate(Ne_RHH_strongsel = 10^log10_popsize * Ne_N_RHH_strongsel) %>%
  mutate(Ne_strongBGS = 10^log10_popsize * Ne_N_HK95_strongBGS) %>%
  mutate(Ne_RHH_BGS_strongsel_strongBGS = 10^log10_popsize * Ne_N_RHH_BGS_strongsel_strongBGS) %>%
  mutate(Ne_RHH_BGS = 10^log10_popsize * Ne_N_RHH_BGS) %>%
  mutate(pi_RHH_BGS_1em9=four_alleles(Ne_RHH_BGS, 1e-9)) %>%
  mutate(pi_RHH_BGS_1em8=four_alleles(Ne_RHH_BGS, 1e-8)) %>%

  mutate(pi_RHH_1em9=four_alleles(10^log10_popsize * Ne_N_RHH, 1e-9)) %>%
  mutate(pi_RHH_1em8=four_alleles(10^log10_popsize * Ne_N_RHH, 1e-8)) %>%

  mutate(pi_BGS_1em9=four_alleles(10^log10_popsize * Ne_N_HK95, 1e-9)) %>%
  mutate(pi_BGS_1em8=four_alleles(10^log10_popsize * Ne_N_HK95, 1e-8)) %>%

  mutate(pi_strongBGS_1em9=four_alleles(Ne_strongBGS, 1e-9)) %>%
  mutate(pi_strongBGS_1em8=four_alleles(Ne_strongBGS, 1e-8)) %>%
   mutate(pi_RHH_1em9_strongsel=four_alleles(Ne_RHH_strongsel, 1e-9)) %>%
  mutate(pi_RHH_1em8_strongsel=four_alleles(Ne_RHH_strongsel, 1e-8)) %>%
 
  mutate(pi_RHH_BGS_1em9_strongsel=four_alleles(Ne_RHH_BGS_strongsel, 1e-9)) %>%
  mutate(pi_RHH_BGS_1em8_strongsel=four_alleles(Ne_RHH_BGS_strongsel, 1e-8)) %>%
  # 4 * U and 10 * J22
  mutate(pi_RHH_BGS_1em9_strongsel_strongBGS=
         four_alleles(Ne_RHH_BGS_strongsel_strongBGS, 1e-9)) %>%
  mutate(pi_RHH_BGS_1em8_strongsel_strongBGS=
         four_alleles(Ne_RHH_BGS_strongsel_strongBGS, 1e-8)) %>%
  filter(phylum %in% unique(da_ml$phylum)) %>%
                      # no point in having one nematode, C. elegans
                      filter(phylum != 'Nematoda')


dml %>% filter(species == "Drosophila melanogaster") %>% 
  select(Ne_N_HK95, Ne_N_RHH_BGS, pi_RHH_BGS_1em8, pi_RHH_BGS_1em9, 
         Ne_RHH_BGS, diversity)  %>%
  mutate_at(vars(Ne_N_HK95, Ne_N_RHH_BGS), function(x) 1-x) %>%
  mutate(Ne_RHH_BGS=signif(Ne_RHH_BGS, 2))

Nc <- dml$log10_popsize
plot(Nc, log10(dml$pi_RHH_BGS_1em9), xlab='', 
     type='n', axes=FALSE,
     ylab='', ylim=c(-3.25, 0),
     xlim=c(4, 15),
     col=phyla_cols[da_ml$phylum], pch=1)
logN <- seq(0, max(Nc), length.out=100)

y1 <- (log10(4e-8) + logN - log10(1 + 4/3 * 4e-8 * 10^logN))
y2 <- (log10(4e-9) + logN - log10(1 + 4/3 * 4e-9 * 10^logN))
polygon(c(logN, rev(logN)), c(y1, rev(y2)), border=NA, 
        col=scales::alpha('gray88', 0.4))
lines(x <- seq(10, 15), rep(log10(3/4), length(x)), col='gray88', 
      lty=2, lwd=0.1, lend=1)

convex_interval <- function(x, y1, y2, n=100, df=10, 
                            col='gray88', alpha=0.2) { 
  h <- alphahull:::ashape(c(x, rev(x)), c(y1, rev(y2)), alpha=1)
  dx <- as.data.frame(h$edges)
  x <- c(dx$x1, dx$x2) 
  sg <- c(dx$y1, dx$y2) 
  # upper line
  uy <- sg[sg %in% y1]
  ux <- x[sg %in% y1]
  uidx <- order(ux)
  uy <- uy[uidx]
  ux <- ux[uidx]

  # lower line
  ly <- sg[sg %in% y2]
  lx <- x[sg %in% y2]
  lidx <- order(lx)
  ly <- ly[lidx]
  lx <- lx[lidx]

  # plot(h)
  # lines(ux, uy, col='red')
  # lines(lx, ly, col='green')
  yl <- smooth.spline(lx, ly, df=df)
  yu <- smooth.spline(ux, uy, df=df)
  col <- scales::alpha(col, alpha)
  polygon(c(yl$x, rev(yu$x)), c(yl$y, rev(yu$y)), border=NA, col=col)
}

idx <- order(dml$log10_popsize)
x <- dml$log10_popsize[idx]
# y1 <- log10(dml$pi_RHH_BGS_1em8[idx])
# y2 <- log10(dml$pi_RHH_BGS_1em9[idx])
# convex_interval(x, y1, y2, col='cornflowerblue')
# xl <- 17.8

y1 <- log10(dml$pi_BGS_1em8[idx])
y2 <- log10(dml$pi_BGS_1em9[idx])
convex_interval(x, y1, y2, col='purple', df=6)
# text(6.1, -1, latex2exp::TeX("$\\pi_{BGS}$"), cex=0.6)

y1 <- log10(dml$pi_strongBGS_1em8[idx])
y2 <- log10(dml$pi_strongBGS_1em9[idx])
convex_interval(x, y1, y2, col='red', df=6)
# text(12, -0.35, latex2exp::TeX("$\\pi_{strong BGS}$"), cex=0.6)

y1 <- log10(dml$pi_RHH_1em8_strongsel[idx])
y2 <- log10(dml$pi_RHH_1em9_strongsel[idx])
convex_interval(x, y1, y2, col='green', df=6)
# text(13, -3.1, latex2exp::TeX("$\\pi_{strong HH}$"), cex=0.6)

y1 <- log10(dml$pi_RHH_1em8[idx])
y2 <- log10(dml$pi_RHH_1em9[idx])
convex_interval(x, y1, y2, col='yellow2', df=6)
# text(11, -0.8, latex2exp::TeX("$\\pi_{HH}$"), cex=0.6)



# text(xl, -0.3, latex2exp:::TeX("$\\mu = 10^{-8}$"), 
     # cex=0.4, col='cornflowerblue', xpd=TRUE)
# text(xl, -0.9, latex2exp:::TeX("$\\mu = 10^{-9}$"), 
     # cex=0.4, col='cornflowerblue', xpd=TRUE)

# y1 <- log10(dml$pi_RHH_BGS_1em8_strongsel[idx])
# y2 <- log10(dml$pi_RHH_BGS_1em9_strongsel[idx])
# convex_interval(x, y1, y2, col='green')
# text(13, -3.1, latex2exp::TeX("$\\pi_{BGS + strong HH}$"), cex=0.6)

# with(dml, segments(log10_popsize, log10(pi_RHH_BGS_1em8), 
#                    log10_popsize, log10(pi_RHH_BGS_1em9),
#                    lwd=0.8,
#                    lend=1,
#                    col=phyla_cols[phylum], pch=1))

# dml <- dml %>%  rowwise() %>%
  # mutate(y = mean(log10(pi_RHH_BGS_1em8), log10(pi_RHH_BGS_1em9)))
# l <- loess(y ~ log10_popsize, dml, span=0.2)
z <- seq(4, 15, length.out=100)
# lines(z, predict(l, data.frame(log10_popsize=z)))
da_dps2 <- da_dps %>% filter(phylum %in% unique(da_ml$phylum)) %>%
                      # no point in having one nematode, C. elegans
                      filter(phylum != 'Nematoda')

# OLS line
phi <- mean(rstan:::extract(dps_fit, pars='phi')$phi)
beta <- mean(rstan:::extract(dps_fit, pars='beta')$beta)
# lines(z, phi + beta*z, lty=2, col='cornflowerblue')
lines(z, coef(dps_ols)[1] + z*coef(dps_ols)[2], lty=2, col='gray52')
# abline(dps_ols, lty=2, col='gray52')


points(da_dps2$log10_popsize, da_dps2$log10_diversity, pch=19, 
       cex=0.8,
       col=phyla_cols[da_dps2$phylum])

title(xlab="approximate population size", line=1.8)
mtext("   diversity (differences per bp)", line=1.9, side=2, xpd=TRUE)


axis(1, seq(4, 14, 2), line=0.3,
     padj = -0.9,
     cex.axis=0.8,
     tck=-0.02,
     labels=latex2exp::TeX(sprintf("$10^{%d}$", seq(4, 14, 2))))
axis(2, seq(0, -3), 
     las=1, 
     cex.axis=0.8,
     tck=-0.02, hadj=0.65,
     labels=latex2exp::TeX(sprintf("$10^{%d}$", seq(0, -3))),
     line=-0.5)
mtext("A", 3, at=4.2, line=-0, cex=1.4, font=2)

legend(4, -0.01, 
       latex2exp::TeX(c("$\\pi_{BGS}$", "$\\pi_{strongBGS}$", 
                        "$\\pi_{HH}$", "$\\pi_{strongHH}$")), 
       fill=alpha(c('purple', 'red', 'yellow2', 'green'), 0.4),
       bty='n', border=0, cex=0.5, ncol=1, xpd=TRUE)

## subfigure B

Nc <- dml$log10_popsize
plot(Nc, log10(dml$pi_RHH_BGS_1em9), xlab='', 
     type='n', axes=FALSE,
     ylab='', ylim=c(-3.25, 0),
     xlim=c(4, 15),
     col=phyla_cols[da_ml$phylum], pch=1)
logN <- seq(0, max(Nc), length.out=100)

y1 <- (log10(4e-8) + logN - log10(1 + 4/3 * 4e-8 * 10^logN))
y2 <- (log10(4e-9) + logN - log10(1 + 4/3 * 4e-9 * 10^logN))
polygon(c(logN, rev(logN)), c(y1, rev(y2)), border=NA, 
        col=scales::alpha('gray88', 0.4))

lines(x <- seq(10, 15), rep(log10(3/4), length(x)), col='gray88', 
      lty=2, lwd=0.1, lend=1)

idx <- order(dml$log10_popsize)
x <- dml$log10_popsize[idx]
y1 <- log10(dml$pi_RHH_BGS_1em8[idx])
y2 <- log10(dml$pi_RHH_BGS_1em9[idx])
convex_interval(x, y1, y2, col='cornflowerblue')
xl <- 17.8
# text(xl, -0.3, latex2exp:::TeX("$\\mu = 10^{-8}$"), 
     # cex=0.4, col='cornflowerblue', xpd=TRUE)
# text(xl, -0.9, latex2exp:::TeX("$\\mu = 10^{-9}$"), 
     # cex=0.4, col='cornflowerblue', xpd=TRUE)
y1 <- log10(dml$pi_RHH_BGS_1em8_strongsel_strongBGS[idx])
y2 <- log10(dml$pi_RHH_BGS_1em9_strongsel_strongBGS[idx])
convex_interval(x, y1, y2, col='orange')

text(11.5, -0.7, latex2exp::TeX("$\\pi_{BGS+HH}$"), cex=0.6)
text(13.4, -3.05, latex2exp::TeX("$\\pi_{strong BGS + strong HH}$"), cex=0.6)



# with(dml, segments(log10_popsize, log10(pi_RHH_BGS_1em8), 
#                    log10_popsize, log10(pi_RHH_BGS_1em9),
#                    lwd=0.8,
#                    lend=1,
#                    col=phyla_cols[phylum], pch=1))

# dml <- dml %>%  rowwise() %>%
  # mutate(y = mean(log10(pi_RHH_BGS_1em8), log10(pi_RHH_BGS_1em9)))
# l <- loess(y ~ log10_popsize, dml, span=0.2)
z <- seq(4, 15, length.out=100)
# lines(z, predict(l, data.frame(log10_popsize=z)))
da_dps2 <- da_dps %>% filter(phylum %in% unique(da_ml$phylum)) %>%
                      # no point in having one nematode, C. elegans
                      filter(phylum != 'Nematoda')

# OLS line
phi <- mean(rstan:::extract(dps_fit, pars='phi')$phi)
beta <- mean(rstan:::extract(dps_fit, pars='beta')$beta)
# lines(z, phi + beta*z, lty=2, col='cornflowerblue')
lines(z, coef(dps_ols)[1] + z*coef(dps_ols)[2], lty=2, col='gray52')
# abline(dps_ols, lty=2, col='gray52')


points(da_dps2$log10_popsize, da_dps2$log10_diversity, pch=19, 
       cex=0.8,
       col=phyla_cols[da_dps2$phylum])

title(xlab="approximate population size", line=1.8)
mtext("   diversity (differences per bp)", line=1.9, side=2, xpd=TRUE)

axis(1, seq(4, 14, 2), line=0.3,
     padj = -0.9,
     cex.axis=0.8,
     tck=-0.02,
     labels=latex2exp::TeX(sprintf("$10^{%d}$", seq(4, 14, 2))))
axis(2, seq(0, -3), 
     las=1, 
     cex.axis=0.8,
     tck=-0.02, hadj=0.65,
     labels=latex2exp::TeX(sprintf("$10^{%d}$", seq(0, -3))),
     line=-0.5)
mtext("B", 3, at=4.2, line=-0, cex=1.4, font=2)

phyla_cols_s <- phyla_cols[unique(dml$phylum)]
legend(6.5, 0.6, names(phyla_cols_s), fill=phyla_cols_s,
       bty='n', border=0, cex=0.5, ncol=3, xpd=TRUE)

if (pdf) dev.off()
par(opar)
 
