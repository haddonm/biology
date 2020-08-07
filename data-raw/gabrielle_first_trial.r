# If confused contact malcolm@haddon.net.au or 0409 941 891

# run the next 4 lines only once )though it wouldn't matter if it was repeated

if (!require(devtools)){install.packages("devtools")} 
devtools::install_github("https://github.com/haddonm/rutilsMH")
devtools::install_github("https://github.com/haddonm/biology")
devtools::install_github("https://github.com/haddonm/makehtml")


# The run proper----------------------------------------------------------------
library(rutilsMH)    # All three available from www.github.com/haddonm.
library(biology)
library(makehtml)

# obviously you would set your own data and results directories.
datadir <- "C:/Users/User/Dropbox/students/gabby/" # note forward slashes in R
# important step, you need to define resdir
resdir <- paste0(datadir,"results")  # you can alter "results" to anything you want

#read in the site tagging data--------------------------------------------------
gdat <- read.csv(filenametopath(datadir,"helio_recapture.csv"),header=TRUE)
head(gdat,20)   # gabbydat

# records by site---------------------------------------------------------------
sites <- sort(unique(gdat$site))  # first two lines especially important
nsites <- length(sites)
printV(sites)
table(gdat$site)


# Now the work starts----------------------------------------------------------
runname <- "InvLonly"   # makeup your own runname
resfile <- setuphtml(resdir=resdir,runname=runname) # for local website. important
# the individual sites-------------------------------------------------------
starttime <- as.character(Sys.time())
final=TRUE

for (i in 1:nsites) {  
  pick <- which(gdat$site == sites[i])
  ans <- fitIL(gdat[pick,],sitename=sites[i],outliers=TRUE)
  # save resulting ans object to resdir
  filen <- filenametopath(resdir,paste0(sites[i],"_model_fit_",runname,".RData")) 
  if (final) save(ans,file=filen)
  #plot the result
  filen <- filenametopath(resdir,paste0(sites[i],"_summary_plots_",runname,".png"))
  plotprep(width=7,height=6,filename=filen,verbose=FALSE)
  plot(ans,ymin=-0.1,xmin=9,outliers=TRUE,nbreaks=30)  
  caption=paste0("Summary plots for ",sites[i]," including the growth curve fit ",
                 "(red dots are outliers), the residuals of the fit, with ",
                 "outliers removed, then the rate of change of DeltaL with Lt, ",
                 "and finally the density distribution of data vs Lt.")
  addplot(filen,resfile,category=sites[i],caption=caption)
  # detailed plot of result
  filen <- filenametopath(resdir,paste0(sites[i],"_detailed_plot_",runname,".png"))
  plotprep(width=7,height=6,filename=filen,verbose=FALSE)
  minx <- getmin(ans$Lt,mult=1.2)
  plotmodelIL(ans,outliers=TRUE,ymin=-0.1,xmin=minx,maintitle=sites[i])
  projtozero(ans)
  caption=paste0("Detailed model fit for ",sites[i])
  addplot(filen,resfile,category=sites[i],caption=caption)
  # save the summary results
  outmat <- summary(ans,console=FALSE)
  filen <- filenametopath(resdir,paste0(sites[i],"_summary_",runname,".csv"))
  if (final) addtable(outmat,filen,resfile,category=sites[i],
                      caption=paste0("Summary results for ",sites[i]))
  
}

endtime <- as.character(Sys.time())
# setup local website of results----------------------------------------
reportlist <- list(
  runname=runname,
  starttime=starttime,
  endtime=endtime
)

runnotes <- "Gabby's data, Inverse Logistic fit only." # optional, could be ""

#  source(filenametopath(sourcedir,"sourcer.R"))
make_html(replist=reportlist,resdir=resdir,width=500,
          openfile=TRUE,runnotes=runnotes,verbose=FALSE,
          packagename="biology")




# Explore the Dose Response curve-----------------------------------------------

library(rutilsMH)    # All three available from www.github.com/haddonm.
library(biology)
library(makehtml)

# obviously you would set your own data and results directories.
datadir <- "C:/Users/User/Dropbox/students/gabby/" # note forward slashes in R
# important step, you need to define resdir
#read in the site tagging data--------------------------------------------------
gdat <- read.csv(filenametopath(datadir,"helio_recapture.csv"),header=TRUE)
head(gdat,20)   # gabbydat

# records by site---------------------------------------------------------------
sites <- sort(unique(gdat$site))  # first two lines especially important
nsites <- length(sites)
printV(sites)
table(gdat$site)


pick <- which(gdat$site == "Maria")
dat <- droplevels(gdat[pick,])


lens <- seq(1,130,1)
maxL <- 5.0
L50 <- 60.0
invl <- invlog(c(maxL,L50,100),lens)
DR <- doseR(c(maxL,L50,4),lens)


initparDR <- function(Lt,DL) {
  pars <- numeric(4)
  minL <- min(Ltin,na.rm=T)
  maxL <- max(Ltin,na.rm=T)
  extent <- maxL - minL
  lim <- minL+0.2*extent
  pick <- which(Ltin < lim)
  if (length(pick) > 1) {
    pars[1] <- mean(DLin[pick])
    pars[4] <- sd(DLin[pick])
  }
  pars[2] <- minL + 0.5*extent
  pars[3] <- minL + 0.9*extent
} # end of initparDR




#' @param x	a vector containing Lt the names must be exact
#' @param y	a vector containing DL the names must be exact
#' @param siteid a numeric site identifier, default = 0
#' @param outliers default = FALSE, should outliers be identified and then omitted?
#' @param sitename an alphanumeric identifier for the site
#' @param ...	the ellipsis is for any remaining parameters
fitDR <- function(x, y, siteid=0,outliers=FALSE,sitename="",...) {
 
  x <- cbind(dat$Lt,dat$DL)
  x <- x[order(x[,1]),]
  siteid=0
  outliers=FALSE
  sitename=dat$site[1]
  
  
  negLDRC <- function(parsin) { # parsin = pars
    pars <- exp(parsin)
    expDL <- doseR(p=pars,x[,1])
    neglogl <- -sum(dnorm(x[,2],expDL,pars[4],log=T))
    return(neglogl)
  } # constant sd with length
  
  negLDRIL <- function(parsin) { # parsin = pars
    expDL <- doseR(p=parsin,x[,1])
    expSD <- parsin[4]/x[,1]
    neglogl <- -sum(dnorm(x[,2],expDL,expSD,log=T))
    return(neglogl)
  } # inverse decrease in sd with length
  
  negLILC <- function(parsin) {  # 
    pars <- exp(parsin)
    expDL <- invlog(p=pars,x[,1])
    neglogl <- -sum(dnorm(x[,2],expDL,pars[4],log=T))
    return(neglogl)
  } # constant sd with length

  negLILIL <- function(parsin) {  # 
    expDL <- invlog(p=parsin,x[,1])
   # expSD <- invlog(p=c(parsin[4],parsin[3],parsin[3]/0.95),x)
    expSD <- parsin[4]/x[,1]
    neglogl <- -sum(dnorm(x[,2],expDL,expSD,log=T))
    return(neglogl)
  } # inverse decrease in sd with length
  
  negLIL <- function(parsin) {
    expDL <- invlog(p=parsin,x[,1])
    expSD <- invlog(p=c(parsin[4],parsin[3],parsin[3]/0.95),x[,1])
    neglogl <- -sum(dnorm(x[,2],expDL,expSD,log=T))
    return(neglogl)
  } # inverse logistic sd with length
  
  
  x <- cbind(dat$Lt,dat$DL)
  x <- x[order(x[,1]),]
  predLt <- seq(6,19,0.2)  
  # dose response
  pars <- c(2.2,9.5,5,3)
 # pars <- log(pars)
  best <- optim(pars,negLDRIL,method="Nelder-Mead",
                hessian=FALSE,
                control=list(trace=0, maxit=1000))
  # inverse logistic
  pars2 <- c(2.1,9.52,17.5,2.5)
 # pars2 <- log(pars2)
  best2 <- optim(pars2,negLILIL,method="Nelder-Mead",
                hessian=FALSE,
                control=list(trace=0, maxit=1000))

  best$par
  best$value
  best2$par
  best2$value 
  abs(best$value - best2$value)
  
  projzero <- function(funk,pars,first,col=3,...) {
    # funk=doseR; pars=exp(best$par); first=x[1,1]; col=4
    Lt <- seq(0.1,first,0.1)
    y <- funk(pars,Lt)
    lines(Lt,y,lwd=2,col=col)
  }
  
  plotprep(width=7,height=4,newdev = FALSE)
  xmax <- getmax(x[,1])
  plot(x[,1],x[,2],type="p",pch=16,cex=1.0,panel.first=grid,xlim=c(0,xmax),
       xlab="Length at Tagging",ylab="Growth Increment")
  lines(predLt,doseR(best$par,predLt),lwd=2,col=4)
  lines(predLt,invlog(best2$par,predLt),lwd=2,col=2)
  legend("topright",legend=c("InvLog","DoseR"),col=c(2,4),lwd=3,bty="n",cex=1.5)
  projzero(doseR,best$par,x[1,1],col=4)
  projzero(invlog,best2$par,x[1,1],col=2)
  
  

  
  plotprep(width=7,height=4,newdev = FALSE)
  xmax <- getmax(x[,1])
  plot(x[,1],x[,2],type="p",pch=16,cex=1.0,panel.first=grid,xlim=c(0,xmax),
       xlab="Length at Tagging",ylab="Growth Increment")
  lines(predLt,doseR(exp(best$par),predLt),lwd=2,col=4)
  lines(predLt,invlog(exp(best2$par),predLt),lwd=2,col=2)
  legend("topright",legend=c("InvLog","DoseR"),col=c(2,4),lwd=3,bty="n",cex=1.5)
  projzero(doseR,exp(best$par),x[1,1],col=4)
  projzero(invlog,exp(best2$par),x[1,1],col=2)
  
  
  
  
  parsin <- initpars(x,y)
  best <- optim(parsin,negLIL,method="Nelder-Mead",
                hessian=FALSE,
                control=list(trace=0, maxit=1000))
  parsin <- best$par
  mod <- nlm(negLIL,parsin,hessian=T, gradtol = 1e-7)
  parsin <- mod$estimate
  MaxDL <- mod$estimate[1]
  L50 <- mod$estimate[2]
  L95 <- mod$estimate[3]
  MaxSig <- mod$estimate[4]
  xout <- NULL  # will contain the list of outliers if one exists
  yout <- NULL
  L50out <- NULL
  L95out <- NULL
  MaxDLout <- NULL
  MaxSigout <- NULL
  if (outliers) {
    L50out <- L50
    L95out <- L95
    MaxDLout <- MaxDL
    MaxSigout <- MaxSig
    expDL <-  invlog(c(MaxDL,L50,L95),x)
    resids <- abs(y - expDL)
    expSD <- invlog(c(MaxSig,L95,(L95/0.95)),x)
    outers <- resids - 2.576*expSD   #99% confidence limits
    pick <- which(outers > 0)
    if ((length(pick) > 0)==TRUE) {
      xout <- x[pick]
      yout <- y[pick]
      x <- x[-pick]
      y <- y[-pick]
    }
    best <- optim(parsin,negLIL,method="Nelder-Mead",
                  hessian=FALSE,
                  control=list(trace=0, maxit=2000))
    parsin <- best$par
    mod <- nlm(negLIL,parsin,hessian=T, gradtol = 1e-7)
    MaxDL <- mod$estimate[1]
    L50 <- mod$estimate[2]
    L95 <- mod$estimate[3]
    MaxSig <- mod$estimate[4]
  }
  Ltrg <- range(x,na.rm=T)
  xmin <- min(Ltrg[1],50)
  xmax <- max(Ltrg[2],180)
  predLt <- seq(xmin,xmax,1)
  predDL <- invlog(c(MaxDL,L50,L95),predLt)
  expDL <-  invlog(c(MaxDL,L50,L95),x)
  resids <- y - expDL
  Nobs <- length(x)
  ans <- list(mod,MaxDL,L50,L95,MaxSig,predLt,predDL,resids,Nobs,x,y,xout,
              yout,L50out,L95out,MaxDLout,MaxSigout,siteid,sitename)
  names(ans) <- c("model","MaxDL","L50","L95","MaxSig","PredLt",
                  "PredDL","resids","Nobs","Lt","DL","OutLt","OutDL","L50out",
                  "L95out","MaxDLout","MaxSigout","siteid","sitename")
  class(ans) <- "IL"
  return(ans)
}




