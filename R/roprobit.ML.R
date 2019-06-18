#' rank ordered probit with maximum simulated likelihood
#' 
#' This function estimates a rank ordered probit model
#' @param formula rank ~ explanatory variables. lower rank = better
#' @param group.ID variable containing the IDs of the decision making units
#' @param data a data.frame
#' @param na.last =TRUE if missing ranks are interpreted as being among the least preferred
#' @param nCores number of parallel cores to use
#' 
#' @examples put example here
#' 
#' the estimation is under the IIA assumption

roprobit.ML <- function(formula, group.ID, data, na.last=T, nCores=1, initparm=NULL) {
  
  rankvar <- all.vars(formula)[1]
  varnames <- attr(terms(formula), "term.labels")
  nCoef <- length(varnames)
  if (!na.last) data <- data[!is.na(data[[rankvar]]),]
  if (is.null(initparm)) initparm <- rep(0, nCoef)
  
  # generate clean ranks that start at 1
  data$rank.clean <- ave(data[[rankvar]], data[[group.ID]], FUN=function(z) rank(z, na.last=T))
  data$rank.clean[is.na(data[[rankvar]])] <- NA
  data <- data[order(data[[group.ID]],data$rank.clean),]
  group.IDs <- unique(data[[group.ID]])
  Ngroups <- length(group.IDs)
  
  # generate difference matrices that take into account partially ranked alternatives
  M <- list()
  Omega <- list()
  k <- 1
  for (i in group.IDs) {
    ROL.i <- data$rank.clean[data[[group.ID]]==i]
    Nranked.i <- sum(!is.na(ROL.i)) # number of non-NA ranks
    Nchoices.i <- length(ROL.i)
    M_i <- matrix(0, Nchoices.i-1, Nchoices.i)
    for (r in 1:(Nranked.i-1)) {
      M_i[r,r] <- -1
      M_i[r,(r+1)] <- 1
    }
    if (Nranked.i<Nchoices.i) {
      for (r in Nranked.i:(Nchoices.i-1)) {
        M_i[r,Nranked.i] <- -1
        M_i[r,(r+1)] <- 1
      }
    }
    M[[k]] <- M_i
    k <- k+1
  }
  
  # generate design matrix
  data$ones <- 1 # needed to trick model.matrix into ignoring NAs of the depvar
  form <- update(formula, ones~.)
  data.split <- split(data, data[[group.ID]])
  X <- lapply(data.split, FUN=function(z) sparse.model.matrix(form, z))
  
  # helper variables
  MX <- list()
  for (i in 1:Ngroups) MX[[i]] <- M[[i]]%*%X[[i]]
  MXmat <- do.call(rbind, MX)
  XXinv <- solve(t(MXmat)%*%MXmat)
  Proj <- XXinv%*%t(MXmat)
  Vcov <- list()
  for (i in 1:Ngroups) Vcov[[i]] <- M[[i]]%*%t(M[[i]])
  lower <- list()
  for (i in 1:Ngroups) lower[[i]] <- rep(-Inf, dim(M[[i]])[1])
  
  # set up likelihood function
  ll <- function(beta) {
    l_i <- 0
    ll <- 0
    for (i in 1:Ngroups) {
      #Vcov <- Vcov / Vcov[1,1]
      l_i <- mvtnorm::pmvnorm(lower=lower[[i]],
                              upper=-as.vector(MX[[i]]%*%beta),
                              sigma=Vcov[[i]],
                              algorithm=GenzBretz())
      ll <- ll+log(l_i)
    }
    attr(ll, 'error') <- NULL
    attr(ll, 'msg') <- NULL
    ll
  }
  
  # optimize
  res <- optim(initparm, ll, control=list(fnscale=-1, trace=1), method='BFGS', hessian=F)
  
  est <- list()
  est$coef <- res$par
  # est$vcov <- tryCatch({
  #     solve(res$hessian)
  #   }, error = function(e) {
  #     write('Hessian matrix not invertible ...', stderr())
  #   }, finally = {
  #     matrix(NA, nCoef, nCoef)
  #   })
  est$optim.res <- res
  class(est) <- "roprobit"
  
  return(est)
}