



library(biology)
library(hplot)
library(codeutils)


x <- 25:175


Wta <- 0.0000562 
Wtb <- 3.161963
	

alta <- 0.00005531
altb <- 3.219760

y <- Wta * x ^ Wtb
alty <- alta * x ^ altb

plotprep(width=9, height=5,newdev=FALSE)
plot1(x,y,lwd=2)
lines(x,alty,lwd=3,col=2)



data(tasab)
table(tasab$site,tasab$sex)
sites <- sort(unique(tasab$site))
nsite <- length(sites)
scenes <- c("Male1","Female1","Male2","Female2")
models <- makelist(scenes)
sexes <- c("M","F")
count <- 0
for (i in 1:nsite) { #  i = 1
  for (s in 1:2) {
    count <- count + 1
    picksex <- which((tasab$sex %in% c(sexes[s], "I")) & (tasab$site == sites[i]))
    models[[count]] <- fitmaturity(tasab[picksex,],
                                   length="length",mature="mature",lower=50,upper=160)
  }
}
str1(models)
str1(models[["Male1"]])


plotprep(width=10, height=8)
parset(plots=c(2,2))
for (i in 1:length(models)) {
  plotmaturity(models[[i]],label=scenes[i],col=2,xmin=0,xmax=0,CI=FALSE,
                           setpar=FALSE,lwd=2)  
}








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



