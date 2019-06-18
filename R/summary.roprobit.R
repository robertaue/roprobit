#' print summary stats
#' 
 summary.roprobit <- function(x) {
   
   su <- data.frame(coef=x$coef, 
                    se=sqrt(diag(x$vcov)))
   su$t <- su$coef/su$se
   su$p <- 2*(1 - pnorm(abs(su$t)))
   su[['2.5%']]  <- x$beta.ci95[1,]
   su[['97.5%']] <- x$beta.ci95[2,]
   
   print(su)
 }