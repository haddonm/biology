library(rutilsMH)
library(maturity)
library(r4maps)


filen <- "C:/Users/User/Dropbox/rcode/maturity/data-raw/maturity_blacklip.csv"

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
pick <- which(ab1$shelllength < 40)
if (length(pick) > 0) { ab1 <- ab1[-pick,] }

take <- c(3,8,10,11,1,7,4,5)
abmat <- ab1[,take]
abmat[,9] <- NULL

pick <- which((abmat[,3] == "M") | (abmat[,3] == "F"))
abmat[pick,9] <- 1
abmat[-pick,9] <- 0
colnames(abmat) <- c("site","year","sex","length","block","month","lat","long","mature")
sites <- sort(unique(abmat$sit_id))   # change here
nsites <- length(sites)
blkyr <- table(abmat$block,abmat$year)

blkyr

blk=c(13)
yr = c(2000,2001)

pick <- which((abmat$block %in% blk) & (abmat$year %in% yr) &
                (abmat$length >= 70) & (abmat$length <= 180))# &
# (abmat$site %in% c(122,155,172,259)))
# (abmat$sit_id %in% c(2,43,258,259,288)))
#   (abmat$sit_id %in% c(2,43,259)))
tmp <- abmat[pick,]
table(tmp$month)

b11 <- droplevels(abmat[pick,c(1,3,4,6,7,8,9)])
table(b11$mature,b11$site)

dist <- getdist(b11)
index <- order(dist[1,])
tmp <- dist[,index]
disto <- tmp[index,]
disto



pick <- which((abmat$block == blk) & (abmat$year %in% yr) &
                (abmat$length >= 50) & (abmat$length <= 180) &
                (abmat$site %in% c(171)))
b11 <- droplevels(abmat[pick,c(1,3,4,6,7,8,9)])
table(b11$mature,b11$site)


#colnames(b11) <- c("site","sex","length","month","lat","lon","mature")
head(b11)
# pick <- which((b11$mature == 1) & b11$length %in% c(118,120))
# b11 <- b11[-pick,]

table(b11$mature,b11$site)
sitenames <- as.numeric(names(table(b11$site)))
sitenames
nsite <- length(sitenames)
table(b11$lat,b11$site)
table(b11$lon,b11$site)
table(b11$month)

nsite=12
getdist(b11)

plotprep(width=7,height=4,newdev=FALSE)
  propm <- as.matrix(tapply(b11$mature,list(b11$length,b11$month),mean,na.rm=TRUE))
  lens <- as.numeric(rownames(propm)) 
  plot1(lens,propm[,1],type="p",cex=1.0,xlabel="Length mm",
        ylabel="Proportion Mature")
if (nsite > 1) for (i in 2:nsite) points(lens,propm[,i],cex=1.0,col=i)


# Simplify
lower=70
upper=160
ab$month <- as.factor(ab$month)
pickL <- which((b11$length >= lower) & (b11$length <= upper))
ab <- droplevels(b11[pickL,])
propm <- tapply(ab$mature,list(ab$length,ab$month),mean,na.rm=TRUE)
if (nsite > 1) {
  lens <- as.numeric(rownames(propm)) 
} else {
  lens <- as.numeric(names(propm)) 
}
# All Sites
ab$site <- as.factor(ab$site) # make site a factor and not a variable
ab$month <- as.factor(ab$month)
smodel <- glm(mature ~ month + length ,family=binomial, data=ab)
outs <- binglm(smodel)
# No site
model <- glm(mature ~ length ,family=binomial, data=ab)
outm <- binglm(model)

plotprep(width=7,height=4,newdev=FALSE)
plot1(lens,propm[,1],type="p",cex=1.0,xlabel="Length mm",
      ylabel="Proportion Mature")
abline(h=c(0.25,0.5,0.75),col=1,lty=3)
for (i in 2:nsite) points(lens,propm[,i],cex=1.0,col=i)

L <- seq(min(ab$length),max(ab$length),1)
pars <- coef(smodel)
term <- exp(pars[1]+ pars[nsite+1] * L)
ys1 <- term/(1+term)   # for site 1
lines(L,ys1,type="l",lwd=2,col=1)  
for (i in 2:nsite) {
  term <- exp((pars[1]+pars[i]) + pars[nsite+1] * L)
  ys2 <- term/(1+term)   # for site 2 with an extra parameter
  lines(L,ys2,type="l",lwd=2,col=i)
}
label <- rownames(outs$coefficients)
first <- paste0("site",sitenames[1])
legend("topleft",legend=c(first,label[2:nsite]),col=c(1:nsite),lwd=3,bty="n")



# No site plot
lower=80
upper=160
pickL <- which((b11$length >= lower) & (b11$length <= upper))
ab <- droplevels(b11[pickL,])
propm <- tapply(ab$mature,ab$length,mean,na.rm=TRUE)
if (nsite > 1) {
  lens <- as.numeric(rownames(propm)) 
} else {
  lens <- as.numeric(names(propm)) 
}
# All Sites
ab$site <- as.factor(ab$site) # make site a factor and not a variable
# No site
model <- glm(mature ~ length ,family=binomial, data=ab)
outm <- binglm(model)

plotprep(width=7,height=4,newdev=FALSE)
plot1(lens,propm,type="p",cex=1.0,xlabel="Length mm",
      ylabel="Proportion Mature")
abline(h=c(0.25,0.5,0.75),col=1,lty=3)
pars <- coef(model)
term <- exp(pars[1]+ pars[2] * L)
ys1 <- term/(1+term)   # for site 1
lines(L,ys1,type="l",lwd=2,col=1)  


