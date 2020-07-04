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


sites <- sort(unique(ab1$site))   # change here
nsites <- length(sites)
sitemth <- table(ab1$month,ab1$site)

sitemth

sit <- "Gardens"

pick <- which((ab1$site == sit) & (ab1$outlier == 0))
sitdat <- droplevels(ab1[pick,])

mthmat <- table(sitdat$month,sitdat$mature)
mthmat
mths <- as.numeric(rownames(mthmat))
mthlab <- c("Jan","Apr","Aug","Nov")


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

