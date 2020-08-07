library(rutilsMH)
library(biology)
#library(r4maps)


# get data and add mature/immature----------------------------------------------
filen <- "C:/Users/User/Dropbox/rcode2/biology/data-raw/maturity_data.csv"

ab <- read.csv(filen,header=TRUE)
dim(ab)
names(ab) <- tolower(names(ab))
names(ab)
ab$sex <- toupper(ab$sex)
dim(ab)
validsex <- c("F","I","M")
pick <- which(ab$sex %in% validsex)
ab1 <- if(length(pick) > 0) { ab[pick,] }
dim(ab1)
ab1$mature <- NA

pick <- which((ab1[,"sex"] == "M") | (ab1[,"sex"] == "F"))
ab1$mature[pick] <- 1
ab1$mature[-pick] <- 0
properties(ab1)








sites <- sort(unique(ab1$site))   # change here
nsites <- length(sites)
sitemth <- table(ab1$month,ab1$site)

sitemth

sit <- c("Gardens","George III Rock")

pick <- which((ab1$site %in% sit) & (ab1$outlier == 0))
sitdat <- droplevels(ab1[pick,])

mthmat <- table(sitdat$season,sitdat$mature,sitdat$site)
mthmat

tapply(sitdat$meatwt,list(sitdat$site,sitdat$season),mean,na.rm=TRUE)


# Fit and plot maturity curves to a single site---------------------------------
plotprep(width=7,height=6,newdev=FALSE)
parset(plots=c(2,2))
for (i in 1:4) { # i=2
  pick <- which(sitdat$month == mths[i])
  mthdat <- droplevels(sitdat[pick,])
  byl <- as.matrix(table(mthdat$length,mthdat$mature))
  lens <- as.numeric(rownames(byl))
  propm <- byl[,2]/rowSums(byl)
  ylabel <- paste0("Proportion Mature _ ",mthlab[i])
  plot(lens,propm,type="p",pch=16,cex=1,col=1,xlab="Length mm",
       ylab=ylabel)
  
 # sitdat$month <- as.factor(sitdat$month)
  smodel <- glm(mature ~ length ,family=binomial, data=mthdat)
  outs <- binglm(smodel)
  L <- seq(min(lens),max(lens),1)
  pars <- coef(smodel)
  term <- exp(pars[1]+ pars[2] * L)
  ys1 <- term/(1+term)   # for site 1
  lines(L,ys1,type="l",lwd=2,col=2)  
}


# All Sites
sit <- "George III Rock"
sit <- "Gardens"

pick <- which((ab1$site == sit) & (ab1$outlier == 0))# & (ab1$month %in% c(1,8)))
sitdat <- droplevels(ab1[pick,])

sitdat$month <- as.factor(sitdat$month)
smodel <- glm(mature ~ month + length ,family=binomial, data=sitdat)
outs <- binglm(smodel)
# calculate the proportion mature by length
byl <- as.matrix(table(sitdat$length,sitdat$mature))
lens <- as.numeric(rownames(byl))
propm <- byl[,2]/rowSums(byl)
mths <- unique(sitdat$month)
nmth <- length(mths)

plotprep(width=7,height=4,newdev=FALSE)
plot(lens,propm,type="p",cex=1.0,xlab="Length mm",
      ylab=paste0("Proportion Mature   ",sit))
abline(h=c(0.25,0.5,0.75),col=1,lty=3)
L <- seq(min(lens),max(lens),1)
pars <- coef(smodel)
term <- exp(pars[1]+ pars[nmth+1] * L)
ys1 <- term/(1+term)   # for site 1
lines(L,ys1,type="l",lwd=2,col=1)  
for (i in 2:nmth) {
  term <- exp((pars[1]+pars[i]) + pars[nmth+1] * L)
  ys2 <- term/(1+term)   # for site 2 with an extra parameter
  lines(L,ys2,type="l",lwd=2,col=i)
}
label <- rownames(outs$out$coefficients)
first <- paste0("month",mths[1])
legend("topleft",legend=c(first,label[2:nmth]),col=c(1:nmth),lwd=3,bty="n",
       cex=1.25)



# Compare sites in January-----------------------------------------------------

pick <- which((ab1$outlier == 0) & (ab1$month == 1))
alldat <- droplevels(ab1[pick,])
unique(alldat$site)

alldat$month <- as.factor(alldat$month)
alldat$site <- as.factor((alldat$site))
amodel <- glm(mature ~ site + length ,family=binomial, data=alldat)
outs <- binglm(amodel)
byl <- as.matrix(table(alldat$length,alldat$mature))
lens <- as.numeric(rownames(byl))
propm <- byl[,2]/rowSums(byl)
nsite <- length(unique(alldat$site))

plotprep(width=7,height=4,newdev=FALSE)
plot(lens,propm,type="p",cex=1.0,xlab="Length mm",
     ylab="Proportion Mature")
abline(h=c(0.25,0.5,0.75),col=1,lty=3)
L <- seq(min(lens),max(lens),1)
pars <- coef(amodel)
term <- exp(pars[1]+ pars[nsite+1] * L)
ys1 <- term/(1+term)   # for site 1
lines(L,ys1,type="l",lwd=2,col=1)  
for (i in 2:nsite) {
  term <- exp((pars[1]+pars[i]) + pars[nsite+1] * L)
  ys2 <- term/(1+term)   # for site 2 with an extra parameter
  lines(L,ys2,type="l",lwd=2,col=i)
}
label <- rownames(outs$coefficients)
first <- "siteGardens"
legend("topleft",legend=c(first,label[2:nsite]),col=c(1:nsite),lwd=3,bty="n",
       cex=1.25)




# pairwise compare sites in January-----------------------------------------------------

pick <- which((ab1$outlier == 0) & (ab1$month == 1) & (ab1$site == "George III Rock"))
alldat <- droplevels(ab1[pick,])
unique(alldat$site)

alldat$month <- as.factor(alldat$month)
alldat$site <- as.factor((alldat$site))
rmodel <- glm(mature ~ length ,family=binomial, data=alldat)
outr <- binglm(rmodel)

pick <- which((ab1$outlier == 0) & (ab1$month == 1) & (ab1$site == "Gardens"))
alldat <- droplevels(ab1[pick,])
unique(alldat$site)

alldat$month <- as.factor(alldat$month)
alldat$site <- as.factor((alldat$site))
gmodel <- glm(mature ~ length ,family=binomial, data=alldat)
outg <- binglm(gmodel)




byl <- as.matrix(table(alldat$length,alldat$mature))
lens <- as.numeric(rownames(byl))
propm <- byl[,2]/rowSums(byl)



# relative condition through year---------------------------------------
sit <- c("Gardens","George III Rock")

pick <- which((ab1$site %in% sit) & (ab1$outlier == 0))
sitdat <- droplevels(ab1[pick,])

mthmat <- table(sitdat$season,sitdat$mature,sitdat$site)
mthmat

tapply(sitdat$meatwt,list(sitdat$site,sitdat$season),mean,na.rm=TRUE)



sites <- sort(unique(sitdat$site))
plotprep(width=7,height=5,newdev=FALSE)
sit <- "Gardens"
seas <- "Autumn"
invar <- "shellwt"
indep <- "length"
pickM <- which(sitdat$season == seas)
allsite <- droplevels(sitdat[pickM,])
xrnge <- range(allsite$length)
pickS <- which(allsite$site == sites[1])
plot(allsite[pickS,indep],allsite[pickS,invar],type="p",pch=16,cex=1.0,
     col=1,xlim=xrnge,panel.first=grid())
points(allsite[-pickS,indep],allsite[-pickS,invar],pch=16,col=2,cex=1.0)


# fix-up outliers in weights --------------------------------------------
pickO <- which(ab1$outlier == 0)
ab2 <- ab1[pickO,]

pick <- which(ab2$length < 30)
ab2[pick,"length"] <- 91; ab2[pick,"outlier"] <- 1
pick <- which(ab2$meatwt > 400)
ab2[pick,"meatwt"] <- 62.6; ab2[pick,"outlier"] <- 1
pick <- which((ab2$length > 145) & (ab2$totwt < 200))
ab2[pick,"outlier"] <- 1
pick <- which((ab2$length < 50) & (ab2$totwt > 25))
ab2[pick,"outlier"] <- 1
pick <- which((ab2$length > 100) & (ab2$shellwt < 25))
ab2[pick,"outlier"] <- 1
pick <- which((ab2$length > 130) & (ab2$meatwt < 50))
ab2[pick,"outlier"] <- 1
pick <- which((ab2$length < 105) & (ab2$totwt > 200))
ab2[pick,"outlier"] <- 1
pick <- which((ab2$length >=80) & (ab2$length < 85) & (ab2$totwt >= 100))
ab2[pick,]
ab2[pick,"outlier"] <- 1



pickO <- which(ab2$outlier == 0)
ab2 <- ab2[pickO,]
xrnge <- range(ab2$length)
plotprep(width=7,height=7,newdev=FALSE)
parset(plots=c(2,2))
invar="totwt"
plot(ab2[,indep],ab2[,invar],type="p",pch=16,cex=1.0,
     col=1,xlim=xrnge,panel.first=grid(),ylab=invar)
invar="shellwt"
plot(ab2[,indep],ab2[,invar],type="p",pch=16,cex=1.0,
     col=1,xlim=xrnge,panel.first=grid(),ylab=invar)
invar="meatwt"
plot(ab2[,indep],ab2[,invar],type="p",pch=16,cex=1.0,
     col=1,xlim=xrnge,panel.first=grid(),ylab=invar)
invar="viscwt"
plot(ab2[,indep],ab2[,invar],type="p",pch=16,cex=1.0,
     col=1,xlim=xrnge,panel.first=grid(),ylab=invar)





















