
#' @title binglm is a cleaner output from the summary of a binomial GLM
#' 
#' @description binglm should really be an S3 method for binomial GLM
#'     summary descriptions. I will make that conversion later. For now
#'     this outputs the information needed in a tighter form. It also
#'     outputs the summary (invisibly).
#'
#' @param x the binomial glm model
#' @param digits the number of digits printed for the coefficients
#'
#' @return the summary of the model (invisibly)
#' @export
#'
#' @examples
#' \dontrun{
#'  # wait on a data set being included into the package
#' }
binglm <- function(x,digits=6) {
  out <- summary(x)
  print(x$call)
  print(round(out$coefficients,digits))
  cat("\nNull Deviance  ",out$null.deviance, " df ",out$df.null,"\n")
  cat("Resid.Deviance ",out$deviance, " df ", out$df.residual,"\n")
  cat("AIC = ",out$aic,"\n\n")
  return(invisible(out))
} #end of binglm







