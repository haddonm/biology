# If confused contact malcolm@haddon.net.au or 0409 941 891

# run the next 4 lines only once )though it wouldn't matter if it was repeated
# 
# if (!require(devtools)){install.packages("devtools")} 
# devtools::install_github("https://github.com/haddonm/rutilsMH")
# devtools::install_github("https://github.com/haddonm/biology")
# devtools::install_github("https://github.com/haddonm/makehtml")
# 

# Explore the Dose Response curve-----------------------------------------------

library(rutilsMH)    # All three available from www.github.com/haddonm.
library(biology)
library(makehtml)

# obviously you would set your own data and results directories.
datadir <- "C:/Users/User/Dropbox/students/gabby/" # note forward slashes in R

#read in the site tagging data--------------------------------------------------
gdat <- read.csv(filenametopath(datadir,"helio_recapture.csv"),header=TRUE)
head(gdat,20)   # gabbydat

# records by site---------------------------------------------------------------
sites <- sort(unique(gdat$site))  # first two lines especially important
nsites <- length(sites)
printV(sites)
table(gdat$site)


# generalize the model fitting--------------------------------------------------  
# first pick your data and define some input data with column 1=Lt and column 2 = DL
printV(sites)
  pick <- which(gdat$site == "Schouten")
  sitedat <- droplevels(gdat[pick,])
  dat <- sitedat[,c("Lt","DL")]
  dat <- dat[order(dat[,1]),]
  dim(dat)
  
  # An alternative sdfunc:
  
  # could use this in place of sdpow; what happens if you do?
  sdconst <- function(p,predDL) {
    return(p[4] * predDL)
  }

  predLt <- seq(0,17,0.2)  #  for plotting results
  
  
  pars <- c(maxDL=8,L50=7,L95=15,sig=0.6,tau=1)
  invl <- fitgrow(p=pars,grow=invlog,sdfunc=sdpow,dat=dat)
  
  pars <- c(maxDL=4,L50=9.5,c=5,sig=0.6,tau=1)
  dr <- fitgrow(p=pars,grow=doseR,sdfunc=sdpow,dat=dat)
  
  outfit(invl,title="Inverse Logistic",backtran = FALSE) # outfit in rutilsMH
  outfit(dr,title="dose response",backtran = FALSE)
  
  
  dr$estimate
  invl$estimate  
  dr$minimum
  invl$minimum 
  abs(dr$minimum - invl$minimum) # Is this greater than 1.92? 
  
  plotprep(width=7,height=4,newdev = FALSE)
  xmax <- getmax(dat[,1])
  plot(dat[,1],dat[,2],type="p",pch=16,cex=1.0,panel.first=grid(),xlim=c(0,xmax),
       xlab="Length at Tagging",ylab="Growth Increment",ylim=c(0,4.5),)
  predDR <- doseR(dr$estimate,predLt)
  predDRsd <- dr$estimate[4] * predDR ^ dr$estimate[5]
  lines(predLt,predDR,lwd=2,col=4)
  predIL <- invlog(invl$estimate,predLt) 
  predILsd <- invl$estimate[4] * predIL ^ invl$estimate[5]
  lines(predLt,predIL,lwd=2,col=2)
  legend("topright",legend=c("InvLog","DoseR"),col=c(2,4),lwd=3,bty="n",cex=1.5)
  lines(predLt,(predIL + 1.96*predILsd),lwd=1,col=2)
  lines(predLt,(predIL - 1.96*predILsd),lwd=1,col=2)
  lines(predLt,(predDR + 1.96*predDRsd),lwd=1,col=4) 
  lines(predLt,(predDR - 1.96*predDRsd),lwd=1,col=4) 
  
  
  
  
  
 



