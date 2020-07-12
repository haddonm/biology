# If confused contact malcolm@haddon.net.au or 0409 941 891

library(rutilsMH)    # All three available from www.github.com/haddonm.
library(biology)
library(makehtml)

# obviously you would set your own data and results directories.
datadir <- "C:/Users/User/Dropbox/students/gabby/" # note forward slashes in R
resdir <- paste0(datadir,"results")

#  source(paste0(datadir,"candidate.R"))  # some extra functions

#read in all the site tagging data
gdat <- read.csv(filenametopath(datadir,"helio_recapture.csv"),header=TRUE)
head(gdat,20)   # gabbydat

# records by site
sites <- sort(unique(gdat$site))
nsites <- length(sites)
printV(sites)
table(gdat$site)



runname <- "InvLonly"
resfile <- setuphtml(resdir=resdir,runname=runname) # for local website
# plot individual sites-------------------------------------------------------
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

reportlist <- list(
  runname=runname,
  starttime=starttime,
  endtime=endtime
)

runnotes <- "Gabby's data, Inverse Logistic fit only."

#  source(filenametopath(sourcedir,"sourcer.R"))
make_html(replist=reportlist,resdir=resdir,width=500,
          openfile=TRUE,runnotes=runnotes,verbose=FALSE,
          packagename="biology")


