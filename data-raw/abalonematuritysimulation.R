



#                      a         b     Lm50       IQm    n
# Autumn_Gardens -14.31522 0.1315631 108.8088 16.700915 196
# Autumn_George  -20.89287 0.1974590 105.8086 11.127496 158
# Spring_Gardens -31.52024 0.2721554 115.8171  8.073419 176
# Spring_George  -30.19454 0.2638436 114.4411  8.327753 178
# Summer_Gardens -12.80808 0.1113119 115.0648 19.739343 211
# Summer_George  -17.96507 0.1604041 111.9988 13.698056 197
# Winter_Gardens -31.22016 0.2791961 111.8216  7.869825 206
# Winter_George  -20.37526 0.1909390 106.7108 11.507471 190

suppressPackageStartupMessages({
  library(hutils)
  library(hplot)
  library(makehtml)
  library(biology)
  library(knitr)
  library(MASS)
})

param <- c(-18.52845, 0.1681592, 110.184, 13.06634)
lens <- 50:150

matur <- maturity(ina=param[1],inb=param[2],lens)


plotprep(width=7, height=4.5, newdev=FALSE)
parset()
plot1(lens,matur,defpar=FALSE)


plotdat <- ans$model$plotdat

plotprep(width=8, height=4.5, newdev=FALSE)
parset()
plotmaturity(ans$model)


L95 <- maturity(CI[1,1],CI[2,1],plotdat[,"lens"])

U95 <- maturity(CI[1,2],CI[2,2],plotdat[,"lens"])

samdat <- lens

predm <- predict(model,newdata=samdat,interval="response")



# My version
model <- glm(mature~length, data=sams, family=binomial(link="logit"))
newdat <- data.frame(length=(50:160))
preddat <- predict(model, newdata=newdat, se.fit=TRUE)

head(data.frame(preddat))




str1(ans)



plotprep(width=8, height=4.5, newdev=FALSE)
plotcomparison(ans=ans,name1=samdat1[1,"siteseas"],name2=samdat2[1,"siteseas"])

















