#' A rank ordered probit model
#' 
#' This function estimates a rank ordered probit model
#' @param formula rank ~ explanatory variables. lower rank = better
#' @param group.ID variable containing the IDs of the decision making units
#' @param data a data.frame
#' @param na.last =TRUE if missing ranks are interpreted as being among the least preferred
#' @param niter number of data augmentation iterations
#' @param thin thinning parameter
#' @param burnin number of discarded samples (in the thinned series)
#' @param method 'MSL' for maximum simulated likelihood, 'Gibbs' for data augmentation
#' @param nCores number of parallel cores to use (only with method='Gibbs')
#' 
#' @examples put example here
#' 


roprobit <- function(formula, group.ID, data, na.last=T, niter=500, thin=10, burnin=50, method='Gibbs', nCores=1) {
  
  # consistency checks
  stopifnot(method %in% c('MSL','Gibbs', 'Gibbs.R'))
  formula.terms <- terms(formula)
  stopifnot(attr(formula.terms, 'response')==1)
  stopifnot(attr(formula.terms, 'intercept')==0)
  
  # preliminaries
  varnames <- attr(formula.terms, "term.labels")
  regform <- update(formula, Y~.) # regression formula for latent valuations
  rankvar <- all.vars(formula)[1]
  nSamples <- floor(niter/thin)
  stopifnot(burnin<=nSamples)
  
  # generate design matrix and bring data in order
  order.vec <- order(data[[group.ID]], data[[rankvar]], na.last=TRUE)
  
  data <- data[order.vec,]
  IDs <- data[[group.ID]]
  nIDs <- length(unique(IDs))
  ROL.length <- aggregate(x=data[[rankvar]], by=list(group.ID=data[[group.ID]]), FUN=function(z) sum(!is.na(z)))[,2]
  if (na.last) { 
    ChoiceSetLength <- as.vector(table(IDs)) 
    # override default behavious to drop rows where the depvar is missing:
    data$zdheiffuj82j <- 1
    X <- sparse.model.matrix(update(formula, zdheiffuj82j~.), data)
    outdata <- data[,c(group.ID, rankvar)]
  } else { 
    ChoiceSetLength <- ROL.length 
    X <- sparse.model.matrix(formula, data) # will omit rows where rank information is missing
    outdata <- data[!is.na(data[[rankvar]]),c(group.ID, rankvar)]
  }
  nCoef <- dim(X)[2]
  MaxUnranked <- rep(-Inf, nIDs)
  MinRanked <- rep(Inf, nIDs)
  
  # simulate ordered probit model (Gibbs sampling, data augmentation)
  XXinv <- solve(t(X)%*%X)
  

  if (method == 'Gibbs.R') {
    # initialize values
    beta <- matrix(0, nrow=nCoef)
    betavalues <- matrix(NA, nrow=nSamples, ncol=nCoef)
    Y <- Xb <- X %*% beta # initialize
    Proj <- XXinv %*% t(X)
    
    for (iter in 1:niter) {
      # simulate latent variable
      k <- 1
      for (i in 1:nIDs) {
        # store variables to reduce memory access
        Nchoices <- ChoiceSetLength[i]
        Nranked <- ROL.length[i]
        MaxUnranked_i <- ifelse(Nranked<Nchoices, Y[(k+Nranked):(k+Nchoices)], -Inf)
        MinRanked_i <- MinRanked[i]
        upper_bound <- Inf
        lower_bound <- -Inf
        for (r in 1:Nchoices) {
          Xb_k <- Xb[k]
          # generate bounds
          if (r==1) {
            upper_bound <- Inf
            lower_bound <- Y[k+1] - Xb_k
          } else if (r<Nranked) {
            upper_bound <- Y[k-1] - Xb_k
            lower_bound <- Y[k+1] - Xb_k
          } else if (r==Nranked) {
            upper_bound <- Y[k-1] - Xb_k
            lower_bound <- MaxUnranked_i - Xb_k
          } else {
            # cat('.', end='')
            upper_bound <- MinRanked_i - Xb_k
            lower_bound <- -Inf
          }
          # draw truncated error term
          u <- ifelse(lower_bound<upper_bound, rtruncnorm(1, a=lower_bound, b=upper_bound, mean=0, sd=1), 0)
          stopifnot(!is.na(u))
          Y[k] <- Xb_k + u
          # update maximum unranked and minimum ranked utility
          if (r==Nranked) {
            MinRanked[i] <- MinRanked_i <- Y[k]
          } 
          # else if (r>Nranked) {
          #   MaxUnranked[i] <- MaxUnranked_i <- max(MaxUnranked[i], Y[k])
          # }
          
          k <- k+1
        }
      }
      # estimate linear model
      
      beta <- Proj %*% Y
      #beta <- rnorm(1, beta.hat, XXinv) # normal prior for beta
      # data$Y_ <- Y_
      # fit <- lm(regform, data=data)
      # beta.hat <- matrix(fit$coefficients, ncol=1)
      if ( !(iter%%thin) ) betavalues[iter/thin,] <- as.matrix(t(beta))
      Xb <- X %*% beta
    }
    
    outdata$latentvalution <- res$Y
    
  } else if (method == 'Gibbs') {
    res <- roprobit_internal(X=X, XXinv=XXinv, niter=niter, thin=thin, ChoiceSetLength=ChoiceSetLength, ROLLength=ROL.length, nCores=nCores)
    betavalues <- res$betadraws
    outdata$latentvalution <- res$Y
  } else if (method == 'MSL') {
    print('Maximum simulated likelihood method not yet implemented.')
    return(0)
  }
  
  # get parameter estimates
  colnames(betavalues) <- colnames(X)
  beta.hat <- colMeans(betavalues[burnin:nSamples,])
  beta.vcov <- cov(betavalues[burnin:nSamples,])
  
  est <- list(coef=beta.hat, vcov=beta.vcov, betavalues=betavalues, niter=niter, burnin=burnin, thin=thin, method=method, valuations=outdata)
  class(est) <- "roprobit"
  
  return(est)
}