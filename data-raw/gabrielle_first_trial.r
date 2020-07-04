# If confused contact malcolm@haddon.net.au or 0409 941 891

library(rutilsMH)    # All three available from www.github.com/haddonm.
library(invLogistic)
library(makehtml)

# obviously you would set your own data and results directories.
datadir <- "C:/Users/User/Dropbox/students/gabby/" # note forward slashes in R
resdir <- paste0(datadir,"results")

source(paste0(datadir,"candidate.R"))  # some extra functions

resfile <- setuphtml(resdir=resdir,runname="FirstR") # for local website

#read in all the site tagging data
gdat <- read.csv(filenametopath(datadir,"growth_increment.csv"),header=TRUE)
head(gdat,20)   # gabbydat

# records by site
sites <- sort(unique(gdat$Site))
nsites <- length(sites)
printV(sites)
table(gdat$Site)


# Fit inverse logistic to all TAS data
# If outliers=FALSE = default, then outliers are ignored.
all <- fitIL(gdat,sitename="AllTAS",outliers=TRUE)
# don't worry about warnings concerning: NANs produced

addtxtobj <- function (object,filen,resfile,category="any",caption = "") {
  nline <- length(object)
  cat(object[1],"\n",file = filen)
  for (i in 2:nline) cat(object[i],"\n",file = filen,append=TRUE)
  logfilename(filename=filen,resfile=resfile,category=category,
              caption=caption)
} # end of addtxtobj




outall <- summary(all) # left column parameters = without outliers
                       # right column = with outliers included
filen <- filenametopath(resdir,"summary_alltas.txt")
addtxtobj(outall,filen,resfile,category="AllTAS",caption="Summary of All TAS fit.")

plotprep(width=7,height=6,newdev=FALSE)
plot(all,miny=-0.1,minx=0,nbreaks=30)


plotprep(width=7,height=6,newdev=FALSE)
plotmodelIL(all,outliers=TRUE,ymin=-0.1)
projtozero(all)


pick <- which(gdat$Site == "Hazards")
pdat <- droplevels(gdat[pick,])
head(pdat)
nrow(pdat)

ans <- fitIL(pdat,sitename=pdat$Site[1],outliers=TRUE)


#ans <- fitDR(p,pdat)
source(paste0(datadir,"candidate.R"))
summary(ans)

filen <- filenametopath(resdir,paste0(pdat$Site[1],"_summary.csv"))
plotprep(width=7,height=6,newdev=FALSE)
plot(ans,miny=-0.1,minx=0,outliers=TRUE,nbreaks=30)

plotprep(width=7,height=6,newdev=FALSE)
plotmodelIL(ans,outliers=TRUE,ymin=-0.1)
projtozero(ans)
x <- seq(0.1,ans$PredLt[1],0.1)
y <- invlog(ans$model$estimate,x)
lines(x,y,lwd=2,col=3)

dyn <- getdyn(ans,initL=1.4)
dyn

