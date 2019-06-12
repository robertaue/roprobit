#' print summary stats
#' 
 summary.roprobit <- function(x) {
   
   su <- data.frame(coef=x$coef, 
                    se=sqrt(diag(x$vcov)))
   su$t <- su$coef/su$se
   su$p <- 2*(1 - pnorm(abs(su$t)))
   
   print(su)
 }