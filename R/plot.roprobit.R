#' show diagnostic plot of an estimation object of class 'roprobit'
#' 

plot.roprobit <- function(x, which=1:min(4,length(x$coef))) {
  #nCoef <- length(x$coef)
  nCoef <- length(which)
  
  par(mfrow=c(nCoef,2))
  for (k in which) {
    plot(x$betavalues[,k], type='l', main=names(x$coef)[k], xlab = 'MCMC samples', ylab = '')
    abline(h=x$coef[k], col="blue", lty=2)
    abline(v=x$burnin)
    plot(density(x$betavalues[x$burnin:dim(x$betavalues)[1],k]), main=names(x$coef)[k])
    abline(v=x$coef[k], col="blue", lty=2)
  }
  par(mfrow=c(1,1))
  
}
