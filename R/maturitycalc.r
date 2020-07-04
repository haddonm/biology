
#' @title binglm is a cleaner output from the summary of a binomial GLM
#' 
#' @description binglm should really be an S3 method for binomial GLM
#'     summary descriptions. I will make that conversion later. For now
#'     this outputs the information needed in a tighter form. It also
#'     outputs the summary (invisibly).
#'
#' @param x the binomial glm model
#' @param digits the number of digits printed for the coefficients
#' @param console should the result be printed to the console, default=TRUE
#'
#' @return the summary of the model (invisibly)
#' @export
#'
#' @examples
#' \dontrun{
#'  # wait on a data set being included into the package
#' }
binglm <- function(x,digits=6,console=TRUE) { # x=smodel; digits=6;console=TRUE
  out <- summary(x)
  pars <- coef(out)[,"Estimate"]
  npar <- length(pars)
  nmod <- npar-1
  columns <- c("a","b","L50","IQ")
  rows <- names(pars[1:nmod])
  models <- matrix(0,nrow=nmod,ncol=4,dimnames=list(rows,columns))
  models[1,] <- c(pars[1],pars[npar],-pars[1]/pars[npar],2*log(3)*pars[npar])
  if (nmod > 1) {
    for (i in 2:nmod) {
      par1 <- pars[1] + pars[i]
      models[i,] <- c(par1,pars[npar],-par1/pars[npar],2*log(3)*pars[npar]) 
    }
  }
  if (console) {  
    print(x$call)
    print(round(out$coefficients,digits))
    cat("\nNull Deviance  ",out$null.deviance, " df ",out$df.null,"\n")
    cat("Resid.Deviance ",out$deviance, " df ", out$df.residual,"\n")
    cat("AIC = ",out$aic,"\n\n")
    print(round(models,digits))
  }
  return(invisible(list(out=out,models=models)))
} #end of binglm




