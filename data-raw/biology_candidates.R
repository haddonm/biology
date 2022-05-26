


# samb=sam; sitecol="siteseas"
analyseall <- function(samb,sitecol="siteseas") {
  sites <- sort(unique(samb[,sitecol]))
  nsites <- length(sites)  
  columns <- c("n","date","propM")
  ans <- as.data.frame(matrix(0,nrow=nsites,ncol=length(columns),
                              dimnames=list(sites,columns)))
  for (i in 1:nsites) {  # i = 1
    picks <- which(samb[,sitecol] == sites[i])
    sams <- samb[picks,]
    n <- length(picks)
    ans[i,] <- c(n,sams[1,"date"],sams[1,"propM"])
  }
  return(invisible(ans))
}



