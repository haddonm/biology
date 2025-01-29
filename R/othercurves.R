



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

#' @title fitexp fits an exponential model
#' 
#' @description fitexp fits an exponential model, such as a weight-at-length
#'     model. Anything that uses a ydat = a.xdat^b form can be fitted by this 
#'     function.
#'
#' @param xdat the values for the independent variable, eg length
#' @param ydat the dependent variable, eg weight
#'
#' @return a list of the model, x and y, being the fitted curve, a and b 
#'     the fitted coefficients, and n the number of observations, and the raw 
#'     data
#' @export
#'
#' @examples
#' x <- seq(50,200,25)
#' y <- c(15,55,130,275,490,790,1220)
#' ans <- fitexp(x,y)
#' cat("a = ", ans$a, " and b = ",ans$b,"\n")
fitexp <- function(xdat,ydat) { # xdat=wasw[,"length"]; ydat=wasw[,"wholewt"]
  n <- length(xdat)
  if (length(ydat) != n)
    stop("Length and Weight data of different lengths \n")
  model <- lm(log(ydat) ~ log(xdat))
  minx <- round(min(xdat))
  maxx <- round(max(xdat))
  coeff <- coef(model)
  x <- minx:maxx  # note x and y has now changed
  a <- exp(coeff[1])
  b <- coeff[2]
  y <- a * x ^ b
  return(invisible(list(model=model,x=x,y=y,a=a,b=b,n=n,xdat=xdat,ydat=ydat)))
} # end fitWtL


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


#' @title fitsiteweights fits each site's weight-at-length curve from an SAU 
#' 
#' @description fitsiteweights fits an exponential curve to the data from each 
#'     site within an sau, along with an analysis of the sau data combined. The
#'     names of the variables used are defined below, they will be natural-log
#'     transofmred within fitexp so the raw values are required here.
#'
#' @param blk the block or SAU containing at least one site
#' @param dat the dataframe containing at least site, length, and wholewt 
#'     columns
#' @param sau name of the statistical blocks being used, default='block'
#' @param site name of the site numbers used, default='site'
#' @param len name of the length variable used, default='length'
#' @param wt name of the weight variable used, default='wholewt'
#'
#' @return a list of the model and x and y the predicted exponential curve
#' @export
#'
#' @examples
#' print("wait for data sets")
fitsiteweights <- function(blk,dat,sau="block",site="site",len="length",
                           wt="wholewt") { 
  # blk=10; dat=wasw;sau="block";site="site";len="length";wt="wholewt"
  pickB <- which(dat[,sau] == blk)
  datb <- dat[pickB,]
  sites <- sort(unique(datb[,site]))
  nsite <- length(sites) 
  resblk <- vector(mode="list",length=(nsite+1)) 
  names(resblk) <- c(paste0(blk,":",sites),paste0(blk,":combined"))
  for (i in 1:nsite) { # i = 1
    picks <- which(datb[,site] == sites[i])
    nobs <- length(picks)
    ans <- fitexp(datb[picks,len],datb[picks,wt])
    resblk[[i]] <- ans 
  } 
  ans <- fitexp(datb[,len],datb[,wt])
  resblk[[nsite+1]] <- ans 
  return(invisible(resblk))
} # end of fitsiteweights

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

#' @title plotsiteweights plots each sites weight-at-length curve from an SAU 
#' 
#' @description plotsiteweights fits an exponential curve to the data from each 
#'     site within an sau and, optionally plots them, along with an analysis 
#'     and plot of the sau data combined. 
#'
#' @param resblk the output from fitsiteweights containing a list of fitexp
#'     results
#' @param cex.lab The font size of the axis labels, default=0.75 = small
#'
#' @return Nothing but a plot is generated
#' @export
#'
#' @examples
#' print("wait for data sets")
plotsiteweights <- function(resblk,cex.lab=0.75) { 
  # resblk=resblk
  nsite <- length(resblk)-1 
  label <- names(resblk)
  comb <- nsite+1
  ans <- resblk[[comb]]
  minx <- getmin(ans$x,mult=1)
  maxx <- getmax(ans$x,mult=1)
  maxy <- getmax(ans$y)
  parset(plots=pickbound((nsite+2)),outmargin = c(1.0,0,0,0),
         margin=c(0.3,0.45,0.05,0.05),byrow=FALSE,cex.lab=cex.lab)
  for (i in 1:nsite) { # i = 1
    ans <- resblk[[i]]
    plot1(ans$xdat,ans$ydat,type="p",cex=1,pch=16,
          maxy=maxy,defpar=FALSE,xlim=c(minx,maxx),ylab=label[i])
    lines(ans$x,ans$y,lwd=2,col=2)
    mtext(paste0("  a = ",round(ans$a,7)),side=3,line=-1.1,adj=0,cex=1.1)
    mtext(paste0("  b = ",round(ans$b,5)),side=3,line=-2.5,adj=0,cex=1.1)
    mtext(paste0("  nobs = ",ans$n),side=3,line=-3.75,adj=0,cex=1.1)
  } 
  ans <- resblk[[comb]]
  plot1(ans$xdat,ans$ydat,type="p",cex=1,pch=16,
        maxy=maxy,defpar=FALSE,xlim=c(minx,maxx),ylab=label[comb])
  lines(ans$x,ans$y,lwd=2,col=2)
  mtext(paste0("  a = ",round(ans$a,7)),side=3,line=-1.1,adj=0,cex=1.1)
  mtext(paste0("  b = ",round(ans$b,5)),side=3,line=-2.5,adj=0,cex=1.1)
  mtext(paste0("  nobs = ",ans$n),side=3,line=-3.75,adj=0,cex=1.1)
  plot1(ans$x,ans$y,lwd=3,col=2,maxy=maxy,ylab="Compare Curves",defpar=FALSE)
  for (i in 1:nsite) {
    lines(resblk[[i]]$x,resblk[[i]]$y,lwd=1,col=1)
  }
  lines(ans$x,ans$y,lwd=3,col=2)
  mtext("Length mm",side=1,line=-0.1,outer=TRUE,cex=1.1)
} # end of fitsiteweights


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
