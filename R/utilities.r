## utilities.r -- utility functions
library(HDInterval)
library(rstan)

wescols <- c("#75BFC1", "#CD2216", "#E8A252", "#D4DAB2", "#6B7181")


lambda_est <- function(fit) {
  sigma_r <- rstan:::extract(fit, pars='sigma_r')$sigma_r
  sigma_p <- rstan:::extract(fit, pars='sigma_p')$sigma_p
  l <- sigma_p / (sigma_p + sigma_r)
  c(mean=mean(l), hdi(l))
}

# bootstrap samples
bs_samp  <- function(x) {
  x[sample(nrow(x), replace=TRUE), ]
}


hpdi <- function(x, alpha) {
  coda::HPDinterval(coda::as.mcmc(x), prob = 1 - alpha)
}

ci_polygon <- function(x, fit=NULL, y=NULL, type='quantile',
                       par = 'y_rep', alpha=0.05, 
                       color='red', transp=0.4, ...) {
  xr <- order(x)
  if (!is.null(fit)) {
    y <- rstan:::extract(fit, pars=par)[[par]]
  }
  stopifnot(!is.null(y))
  if (type == 'quantile') {
    cis <- apply(y, 2, function(p) 
                 quantile(p, probs=c(alpha/2, 1-alpha/2)))
  } else if (type == 'hpdi') {
    cis <- apply(y, 2, function(p) hpdi(p, alpha))
  } else {
    stop('type must be hpdi or quantile')
  }
  polygon(c(x[xr], rev(x[xr])), c(cis[2, xr], rev(cis[1, xr])), 
          border=0, col=alpha(color, transp), ...)
  invisible(cis)
}

maplength_bodysize_plot <- function(x, default_axes=FALSE) {
  # wraps global data 
  with(x, plot(map_length, log10_size, pch=19, cex=0.7, col='gray48'))
 if (!default_axes) {
    with(x, plot(map_length, log10_size, pch=19, axes=FALSE, cex=0.7,
                         ylab='', xlab='', 
                         xlim=c(0, 50), ylim=c(-4.3, 1.1), 
                         col='gray48'))
    axis(1, seq(0, 50, 10))
    axis(2, seq(-4, 1), seq(-4, 1), labels=latex2exp::TeX(sprintf("$10^{%d}$", seq(-4, 1))))
  }
}


TAX_GROUPS <- c('kingdom', 'superphylum', 'phylum', 
                'subphylum', 'class', 'order', 'family', 
                'genus', 'species')

extract_rank <- function(x, rank = 'Kingdom', default=NULL) {
  if (!length(x) || is.na(x))
    return(NA)
  qry <- rank %in% x$rankname
  if (!any(qry)) {
    if (!is.null(default)) return(default)
    stop('rank not found')
  }
  x$taxonname[x$rankname == rank[which(qry)]]
}

extract_genus <- function(x, genus_only=FALSE) {
  sub('([^_]+) ([^_]+).*', '\\1',  x)
}



extract_species <- function(x, genus_only=FALSE) {
  replace <- '\\1 \\2'
  if (genus_only) {
    replace <- '\\1'
  }
  sub('([^_]+)_([^_]+).*', replace,  x)
}


# some Leffler names have multiple species like Genus species1/species2
# etc. This just grabs the first one (ok for tree purposes)
clean_name <- Vectorize(function(x) strsplit(x, '/', fixed=TRUE)[[1]][1])


# confidence intervals
ci <- function(x, est, alpha=0.05, pivot=TRUE, with_est=TRUE) {
  cis <- quantile(x, probs=c(alpha / 2, 1-alpha/2))
  if (pivot) {
    if (!with_est)
      return(c(2*est - cis[2], 2*est - cis[1]))
    return(c(2*est - cis[2], est, 2*est - cis[1]))
  }
  if (!with_est)
    return(cis)
  return(cis[1], est, cis[2])
}

# bootstrap confidence interval polygon (a shaded line for 
# figures with posterior densities to compare bayes/freq)
# height is polygon height (x positions are cis)
bs_ci_polygon <- function(ci, height, alpha=0.2) {
 polygon(c(ci[1], ci[1], ci[3], ci[3]), c(0, height, height, 0), 
          border=NA, 
          col=alpha('gray42', alpha))
}

# posterior plot, with optional CIs
post_plot <- function(post, alpha=0.05, main='', thresh=0.90,
                      ylab='', cex.main=1, bs_cis=NULL,
                      title_line=1, cex.lab=1.3,
                      col='cornflowerblue',
                      xlim=NULL) {
  dns <- density(post)
  plot(dns, axes=FALSE, type='n',
       xlab='', 
       ylab='', main=main, cex.main=cex.main, xlim=xlim,
       yaxs="i", ylim=c(0, 1.05*max(dns$y)))
  hpdi <- hdi(post, thresh)
  x <- dns$x
  y <- dns$y
  idx <- hpdi[1] <= x & hpdi[2] >= x 
  xx <- x[idx]
  yy <- y[idx]
  yy[1] <- 0
  yy[length(yy)] <- 0
  xseq <- seq(min(xx), max(xx), length.out=length(yy))

  if (!is.null(bs_cis)) {
    bs_ci_polygon(bs_cis, 1.1*max(yy)) 
  }
  polygon(xseq, yy, border=NA, 
          col=alpha(col, 0.4))
  if (!is.null(bs_cis)) {
    abline(v=bs_cis[2], col='gray22', lend=1)
  }
  ave <- mean(post) 
  segments(ave, 0, ave, y[which.min(abs(x - ave))], 
        col=col, lwd=2, lend=1)
  abline(v=0, lty=2, col='gray42', lwd=0.8, lend=1)
  lines(dns) 
  title(ylab=ylab, line=title_line, cex.lab=cex.lab,
        xpd=TRUE)
  axis(1, tck=-0.08, padj=-1.8, cex.axis=0.7)
  # axis(2)
  return(hpdi)
}


# from phytools, modified for log values
add.color.bar.log<-function(leg,cols,title=NULL,lims=c(0,1),digits=1,prompt=TRUE,lwd=4,outline=TRUE,...){
	if(prompt){
		cat("Click where you want to draw the bar\n")
		flush.console()
		x<-unlist(locator(1))
		y<-x[2]
		x<-x[1]
	} else {
		if(hasArg(x)) x<-list(...)$x
		else x<-0
		if(hasArg(y)) y<-list(...)$y
		else y<-0
	}
	if(hasArg(fsize)) fsize<-list(...)$fsize
	else fsize<-1.0
	if(hasArg(subtitle)) subtitle<-list(...)$subtitle
	else subtitle<-NULL
	if(hasArg(direction)) direction<-list(...)$direction
	else direction<-"rightwards"
	if(direction%in%c("rightwards","leftwards")){
		X<-x+cbind(0:(length(cols)-1)/length(cols),1:length(cols)/length(cols))*(leg)
		if(direction=="leftwards"){ 
			X<-X[nrow(X):1,]
			if(!is.null(lims)) lims<-lims[2:1]
		}
		Y<-cbind(rep(y,length(cols)),rep(y,length(cols)))
	} else if(direction%in%c("upwards","downwards")){
		Y<-y+cbind(0:(length(cols)-1)/length(cols),1:length(cols)/length(cols))*(leg)
		if(direction=="downwards"){ 
			X<-X[nrow(X):1,]
			if(!is.null(lims)) lims<-lims[2:1]
		}
		X<-cbind(rep(x,length(cols)),rep(x,length(cols)))
	}
	if(outline) lines(c(X[1,1],X[nrow(X),2]),c(Y[1,1],Y[nrow(Y),2]),lwd=lwd+2,lend=2) 
	for(i in 1:length(cols)) lines(X[i,],Y[i,],col=cols[i],lwd=lwd,lend=2)
	if(direction%in%c("rightwards","leftwards")){
		if(!is.null(lims)) text(x=x,y=y,
			latex2exp:::TeX(sprintf("$10^{%d}$", as.integer(round(lims[1],digits)))),pos=3,cex=fsize)
		if(!is.null(lims)) text(x=x+leg,y=y,
			latex2exp:::TeX(sprintf("$10^{%d}$", as.integer(round(lims[2],digits)))),pos=3,cex=fsize)
		if(is.null(title)) title<-"P(state=1)"
		text(x=(2*x+leg)/2,y=y,title,pos=3,cex=fsize)
		if(is.null(subtitle)) 
			text(x=(2*x+leg)/2,y=y,paste("length=",round(leg,3),sep=""),pos=1,cex=fsize)
		else text(x=(2*x+leg)/2,y=y,subtitle,pos=1,cex=fsize)
	} else if(direction%in%c("upwards","downwards")){
		if(!is.null(lims)) text(x=x,y=y-0.02*diff(par()$usr[3:4]),round(lims[1],digits),
			pos=1,cex=fsize)
		if(!is.null(lims)) text(x=x,y=y+leg+0.02*diff(par()$usr[3:4]),
			round(lims[2],digits),
			pos=3,cex=fsize)
		if(is.null(title)) title<-"P(state=1)"
		text(x=x-0.04*diff(par()$usr[1:2]),y=(2*y+leg)/2,title,
			pos=3,cex=fsize,srt=90)
		if(is.null(subtitle)) 
			text(x=x+0.04*diff(par()$usr[1:2]),y=(2*y+leg)/2,
				paste("length=",round(leg,3),sep=""),pos=1,
				srt=90,cex=fsize)
		else text(x=x+0.04*diff(par()$usr[1:2]),y=(2*y+leg)/2,
			subtitle,pos=1,cex=fsize,srt=90)
	}
}


