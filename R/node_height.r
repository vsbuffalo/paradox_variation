

getTips <- function(phy, x) {
  tips <- seq_along(phy$tip.label)
  phy$tip.label[intersect(getDescendants(phy, x), tips)]
}

node_height <- function(x, phy, abs_pics=TRUE, log_pics=FALSE, 
                        scaled=TRUE, method='rlm') {
  fit_func <- list('lm'=lm, 'rlm'=MASS:::rlm)[[method]]
  x <- setNames(x, phy$tip.label)
  ics <- pic(x, phy, scaled=scaled)
  if (abs_pics) 
    ics <- abs(ics)
  if (log_pics)
    ics <- log(pics)
  phy$node.label <- NULL  # otherwise we can't access descendants
  bt <- branching.times(phy)
  # this just reverses the x-axis, which I don't like
  # bt <- abs(bt - max(bt))
  d <- tibble(bt= bt, pic=ics, nodes=names(bt))
  d <- d %>% mutate(tips = map(nodes, ~ getTips(phy, .)))
  fit <- fit_func(ics ~ bt)

  df <- summary(fit)$df[2]
  tval <- summary(fit)$coefficients[2, 3]
  pval <- dt(abs(tval), df)
 
  list(d=d, fit = fit, abs_pics=abs_pics, log_pics=log_pics,
      method=method, pval=pval)
}


# plot node-height tests
plot_nh <- function(x, title="", ylim=NULL, 
                    xlab="branching time", twotail=FALSE,
                    lowess_line=FALSE,
                    ylab=NULL, alpha=0.05, f=1/4, at=NULL) {
    lab <- "contrasts"
  if (x$abs_pics)
    lab <- "absolute value of contrasts"
  if (x$log_pics)
    lab <- sprintf('log(%s)', lab)
  if (!is.null(ylab))
    lab <- ylab
  plot(pic ~ bt, x$d, axes=FALSE,
       type='n',
       main='',
       ylim=ylim,
       xlab='', ylab='')
  points(pic ~ bt, x$d, pch=21, 
         cex=1.3,
       bg='gray22', col='white',
       lwd=0.2)
  if (x$method == 'rlm') {
    df <- summary(x$fit)$df[2]
    tval <- summary(x$fit)$coefficients[2, 3]
    pval <- dt(abs(tval), df)
    print(pval)
    if (pval <= alpha/(1 + as.integer(twotail)))
      abline(x$fit, col='red', lwd=1.7)
  }  else {
    stop('not implemented')
  }
  # text(100, 0.5, pval)
  title(ylab=lab, 
        line=1.9, cex.lab=1)
  title(xlab=xlab,
        line=2.5, cex.lab=1.2)
  title(main=title, cex.main=1.4, font.main = 1, line=-1)

  if (pval <= alpha/(1 + as.integer(twotail))) {
    if (lowess_line) {
      lines(lowess(x$d$pic ~ x$d$bt, f=f), 
            lwd=1.7, 
            col='cornflowerblue')
    }
  }
  axis(1, padj = -1.8,
     tck=-0.02,
     cex.axis=0.7, line=0.3)
  axis(2, at=at,
     las=1,
     tck=-0.02, hadj=0.5,
     cex.axis=0.7, line=0.3)

}


