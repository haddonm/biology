



#' @title classicalL the classical logistic curve
#'
#' @description classicalL this uses the logistic function:
#'   exp(a+b*L)/(1+exp(a+b*L)), which has the property that the L50 = -a/b. We
#'   use this for growth increments and use a GLM with binomial errors to 
#'   obtain the same equation for maturity curves (taht assume maxDL=1.0).
#'
#' @param p a vector of three parameters, MaxDL, the maximum growth increment,
#'     'a', the intercept of the exponential function if positive and b is
#'     negative then we have an inverse logistic curve, and 'b' the gradient of 
#'     the exponential function. If 'b' is positive and a is negative then a 
#'     positive logistic curve is generated.
#' @param lens a vector of lengths for which the logistic value
#'     will be calculated
#' @return A vector of length(lens) containing the predicted logistic value at
#'     length values (can be used for growth increments, we use a statistical 
#'     method to estimate the same equation for maturity, with MaxDL = 1.0)
#' @export
#'
#' @examples
#' maxDL <- 5.0
#' a <- 6
#' b <- -0.1
#' lens <- seq(1,130,1)
#' growthinc <- classicalL(c(maxDL,a,b),lens)
classicalL <- function(p,lens) { # p=c(6,7,-1.5,1,1); lens=
  ans <- p[1]*exp(p[2]+p[3]*lens)/(1+exp(p[2]+p[3]*lens))
  return(ans)
} # end of classicalL

#' @title doseR the dose response estimate of deltaL
#'
#' @description doseR calculates the dose response estimates of deltaL growth
#'     from input parameters a, b, and c: deltaL <- a/(1 + ((L/b)^c))
#'
#' @param p the parameter list, only the first three value are used here.
#'     Typically p would also contain the sigma value representing the standard
#'     deviation of the additive normal random residuals expected about the
#'     fitted curve.
#' @param lens the initial or starting lengths
#'
#' @return a vector of the same length as nrow(x)
#' @export
#'
#' @examples
#' p=c(1.46,11.0,6.5,0.25) # a, b, c, sigma, the latter not used here
#' # x=pdat
#' print("wait on internal data")
doseR <- function(p,lens) {#
  ans <- p[1]/(1 + ((lens/p[2])^p[3]))
  return(ans)
} # end of doseR

#' @title fitgrow fits a tagging growth model to tagging data
#' 
#' @description fitgrow fits a selected tagging growth model to a set of 
#'     initial lengths and consequent growth increments. It requires the growth 
#'     model function to be input as a function argument as well as another 
#'     function, as a second argument, that is used to describe how the variance
#'     is expected to change with the changing predicted growth increments. The
#'     function uses negLLG to calculate the negative log-likelihoods
#'      
#' @param p the parameters a, b, c, sigma, and tau
#' @param grow the function describing the growth function used
#' @param sdfunc the function describing the sd vs predicted DL
#' @param dat a data.frame containing at least Lt (initial length) and DL (delta
#'     L) in columns 1 and 2
#' @param hessian should the hessian be generated? default=FALSE
#' 
#' @return the fitted model output from the nlm function 
#' @export
#'
#' @examples
#' print("waiting to be developed")     
fitgrow <- function(p,grow,sdfunc,dat,hessian=FALSE) { # p=c(1.46,12.0,6.5,0.5); dat=pdat;Lt="Lt";DL="DL"
  best <- optim(p,negLLG,grow=grow,sdfunc=sdfunc,x=dat,method="Nelder-Mead",
                hessian=hessian,
                control=list(trace=0, maxit=1000))
  mod <- nlm(negLLG,best$par,grow=grow,sdfunc=sdfunc,x=dat,hessian=hessian, 
             gradtol = 1e-7)
  return(mod)
} # end of fitgrow

#' @title negLLG calculates the negative log-likelihood for growth models
#' 
#' @description negLLG is a general function for calculating the negative log-
#'     likelihood using additive normal random errors for tagging based growth 
#'     models defined in the functions 'grow' and 'sdfunc'. The 'grow' 
#'     argument must be a function that defines the expected growth increments 
#'     given two input arguments 'p' and 'Lt', where p is the vector of 
#'     parameters, and 'x' (or whatever) is a vector of initial lengths. The
#'     sdfunc is a function that defines the description of the standard 
#'     deviation around the predicted growth increments. A penalty is imposed on
#'     parameter[2], which is often L50 to avoid negative numbers
#' @param p the vector of parameters, what these are will depend on both grow
#'     and sdfunc 
#' @param grow the function describing the growth function used
#' @param sdfunc the function describing the sd vs predicted DL
#' @param x a data.frame or matrix containing at least Lt (initial length) 
#'     and DL (deltaL) in columns 1 and 2
#' @param pen a penalty to keep p[2] positive. default = TRUE
#'     
#' @return a single number, the negative log-likelihood
#' @export    
#'     
#' @examples 
#'  pars <- c(2.2,9.5,5,0.6,1)
#'  dat <- NULL  # wait on internal dataset
#'  # negLLG(p=pars,grow=invlog,sdfunc=sdpow)
negLLG <- function(p,grow,sdfunc,x,pen=TRUE) { # par=p; x=pdat
  expDL <- grow(p=p,x[,1])
  expSD <- sdfunc(p=p,expDL)
  neglogl <- -sum(dnorm(x[,2],expDL,expSD,log=T))
  if (pen) neglogl <- neglogl + 100*exp(-1000*p[2])
  return(neglogl)
}

#' @title outIL provides a tidy summary of the output from fitIL
#'
#' @param object the object output from fitIL
#' @param console should the output be printed to the console, default=TRUE
#'
#' @return a character vector of 16 statements ready to be printed
#' @export
#'
#' @examples
#' print("waiting to be developed")
outIL <- function(object,console=TRUE) {
  txt <- vector("character",16)
  Ltrge <- range(object$Lt,na.rm=TRUE)
  DLrge <- range(object$DL,na.rm=TRUE)
  outs <- FALSE
  if (length(object$MaxDLout) > 0) { outs <- TRUE }
  txt[1] <- paste0("siteid  : ",object$siteid)
  txt[2] <- paste0("sitename: ",object$sitename)
  if (outs){ txt[3] <- paste0("MaxDL   : ",round(object$MaxDL,digits=4),"   ",
                              round(object$MaxDLout,digits=4))
  } else { txt[3] <- paste0("MaxDL   : ",round(object$MaxDL,digits=4))
  }
  if (outs){ txt[4] <- paste0("L50     : ",round(object$L50,digits=4),"  ",
                              round(object$L50out,digits=4))
  } else { txt[4] <- paste0("L50     : ",round(object$L50,digits=4))
  }
  if (outs){ txt[5] <- paste0("L95     : ",round(object$L95,digits=4)," ",
                              round(object$L95out,digits=4))
  } else { txt[5] <- paste0("L95     : ",round(object$L95,digits=4))
  }
  if (outs){ txt[6] <- paste0("MaxSig  : ",round(object$MaxSig,digits=4),"   ",
                              round(object$MaxSigout,digits=4))
  } else { txt[6] <- paste0("MaxSig  : ",round(object$MaxSig,digits=4))
  }
  txt[7] <- paste0("N       : ",object$Nobs)
  txt[8] <- paste0("Outliers: ",length(object$OutLt))
  txt[9] <- paste0("Range Lt: ",Ltrge[1],"    ",Ltrge[2])
  txt[10] <- paste0("Range DL: ",DLrge[1],"    ",DLrge[2])
  txt[11] <- paste0("-ve LL  : ",round(object$model$minimum,5))
  txt[12] <- paste0("Other Components")
  txt[13] <- paste0("$Lt and $DL are the input data minus any outliers")
  txt[14] <- paste0("$model  contains the nlm fit")
  txt[15] <- paste0("$PredLt and PredDL = fitted line")
  txt[16] <- paste0("$OutLt and $OutDL = outlier values")
  if (console) for (i in 1:16) cat(txt[i],"\n")
  return(txt)
} # end of outIL



#' @title projtozero extends the predicted inverse logistic to zero
#'
#' @description projtozero extends the predicted inverse logistic derived from
#'     fitIL to zero.
#'
#' @param ans the output from fitIL
#'
#' @return nothing but does add a line to a plot from plotmodelIL
#' @export
projtozero <- function(ans) {
  x <- seq(0.1,ans$PredLt[1],0.1)
  y <- invlog(ans$model$estimate,x)
  lines(x,y,lwd=2,col=3)
} #end of projtozero

#' @title sdpow a power relationship used to describe the changing variance
#' 
#' @description sdpow a power relationship p[4]*predDL^p[5] used to describe
#'     how the standard deviation to be used when estimating normal likelihoods
#'     changes relative to the predicted growth increment
#'
#' @param p the vector of parameters, what these are will depend on both grow
#'     and sdfunc 
#' @param predDL the predicted DeltaL derived from some growth function
#'
#' @return a vector of standard deviations the same length as predDL
#' @export
#'
#' @examples
#'  pars <- c(2.2,9.5,5,0.6,1)
#'  dat <- NULL  # wait on internal dataset
#'  # negLLG(p=pars,grow=invlog,sdfunc=sdpow)
sdpow <- function(p,predDL) {
  return(p[4] * predDL ^ p[5])
}
