

library(rutilsMH)
library(biology)

options("show.signif.stars"=FALSE,"stringsAsFactors"=FALSE,
        "max.print"=50000,"width"=240)



data(midg)

ans <- fitIL(midg$Lt,midg$DL,outliers=TRUE,sitename="Middle Ground")

print(ans)

plotprep(width=7,height=6,newdev=FALSE)
plot(ans)

plotprep(width=7,height=5,newdev=FALSE)
plotmodelIL(ans,outliers=TRUE)
projtozero(ans)

dyn <- getdyn(ans,maxage=30)
dyn

summary(ans)


boot <- dobootIL(ans,reps=100)

plotprep(width=7,height=6,newdev=FALSE)
plot(boot,col=4)




