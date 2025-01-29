
#' @title analysegroups identifies tentative groups of sites in a given block
#' 
#' @description analysegroups uses a block and a data.frame of maturity data to
#'     conduct analyses that characterize the properties of all sites within a
#'     block that are present in the data.frame. It fits maturity ogives to each
#'     site, and identifies blocks based on attempts to combine sites in the
#'     sequential order of their latitudes, first low to high and then the 
#'     reverse.  
#'
#' @param blk a statistical block or SAU containing sites of maturity samples 
#' @param samdat a data.frame containing the raw maturity data used in the 
#'     model fitting
#' @param sitecol column name of the variable for the site numbers. The default
#'     column name = 'site'
#'
#' @return a list of nsites, siteproperties, site maturity parameters, the 
#'     fitted models, the tentative groups, and the selected data fro mteh input 
#'     data.frame
#' @export
#'
#' @examples
#' args(analysegroups)  #  ?analysegroups
#' print("wait on internal datasets")
#' # blk=5; samdat=samnw; sitecol="site"
analysegroups <- function(blk,samdat,sitecol="site") {
  pickB <- which(samdat$block == blk)
  samb <- samdat[pickB,]
  sites <- sort(unique(samb[,sitecol]))
  nsites <- length(sites)  
  sprops <- siteprops(samb)
  outsite <- getsiteparams(samb)
  sitepar <- outsite$sitepar  
  if (length(sites) > 1) {
    sitepar <- sitepar[order(sitepar[,"Lm50"]),]
    sitenums <- as.numeric(rownames(sitepar))
    #estimate group structure
    groups1 <- findgroups(sitenums=sitenums,samdat=samb)
    #reverse the ordering  on Lm50
    sitepar2 <- outsite$sitepar
    sitepar2 <- sitepar2[order(sitepar2[,"Lm50"],decreasing=TRUE),]
    sitenums <- as.numeric(rownames(sitepar2))
    nsites <- length(sitenums)
    groups2 <- findgroups(sitenums=sitenums,samdat=samb)
  } else {
    groups1 <- list(sites)
    groups2 <- NULL
  }
  return(invisible(list(nsites=nsites,sprops=sprops,sitepar=sitepar,
                        outsite=outsite,group1=groups1,group2=groups2,
                        samb=samb)))
} # end of analysegroups

#' @title binglm is a cleaner output from the summary of a binomial GLM
#' 
#' @description binglm should really be an S3 method for binomial GLM
#'     summary descriptions. I will make that conversion later. For now
#'     this outputs the information needed in a tighter form. It also
#'     outputs the summary (invisibly).
#'
#' @param x the binomial glm model
#' @param digits the number of digits printed for the coefficients
#' @param console should the result be printed to the console, default=TRUE
#'
#' @return the summary of the model and the parameters (invisibly)
#' @export
#'
#' @examples
#' args(binglm)   #   ?binglm
#'  # wait on a data set being included into the package
binglm <- function(x,digits=6,console=TRUE) { # x=smodel; digits=6;console=TRUE
  out <- summary(x)
  pars <- coef(out)[,"Estimate"]
  npar <- length(pars)
  nmod <- npar-1
  columns <- c("a","b","L50","IQ")
  rows <- names(pars[1:nmod])
  models <- matrix(0,nrow=nmod,ncol=4,dimnames=list(rows,columns))
  models[1,] <- c(pars[1],pars[2],-pars[1]/pars[2],2*log(3)/pars[2])
  if (nmod > 1) {
    for (i in 2:nmod) {
      par1 <- pars[1] + pars[1+i]
      models[i,] <- c(par1,pars[2],-par1/pars[2],2*log(3)/pars[2]) 
    }
  }
  if (console) {  
    print(x$call)
    print(round(out$coefficients,digits))
    cat("\nNull Deviance  ",out$null.deviance, " df ",out$df.null,"\n")
    cat("Resid.Deviance ",out$deviance, " df ", out$df.residual,"\n")
    cat("AIC = ",out$aic,"\n\n")
    print(round(models,digits))
  }
  return(invisible(list(out=out,param=models)))
} #end of binglm

#' @title bootmature generates bootstrap 95p CI for a maturity logistic curve  
#' 
#' @description bootmature takes the data used when fitting a maturity logistic
#'     curve, conducts 1000 bootstrap samples and refits the model to obtain 
#'     the parameters and predicted probabilities of maturity.
#'
#' @param samdat the maturity dataset with a minimum of length and maturity
#'     identified by the arguments 'lens' and 'mature'.
#' @param label a character string for a label to the model fit, default=''
#' @param lens the column name of the length data in the data.frame, default =
#'     'length'
#' @param mature the column name of the maturity data 1 = mature, 0 = immature,
#'     default='mature'
#' @param nb the number of bootstraps to perform, default=1000
#' @param lower the smallest size class for inclusion in the predicted plot 
#'     data, the default = 50
#' @param upper the largest size class for inclusion in the predicted plot 
#'     data, the default = 160 
#'     
#' @seealso {
#'   \link{fitmaturity}, \link{sample}
#' }
#'
#' @return a list of bootpar, bootout, and bootCI
#' @export
#'
#' @examples
#' print("wait on data sets")
#' # samdat=dat; label="";lens="length";mature="mature";nb=1000;lower=50;upper=160
bootmature <- function(samdat,label="",lens="length",mature="mature",nb=1000,
                       lower=50, upper=160) {
  dat <- samdat[,c(lens,mature)] 
  nr <- nrow(dat)
  index <- 1:nr
  predL <- lower:upper
  nL <- length(predL)
  bootout <- matrix(0,nrow=nb,ncol=nL,dimnames=list(1:nb,predL))
  columns <- c("a","b","L50","IQ")
  bootpar <- matrix(0,nrow=nb,ncol=4,dimnames=list(1:nb,columns))
  for (i in 1:nb) {
    samp <- sample(index,nr,replace=TRUE)
    boot <- dat[samp,]
    ans2 <- fitmaturity(boot,length=lens,mature=mature,lower=lower,upper=upper)
    bootout[i,] <- ans2$plotdat[,"predm"]
    bootpar[i,] <- ans2$param
  }
  bootCI <- apply(bootout,2,quants)
  return(invisible(list(bootpar=bootpar,bootout=bootout,bootCI=bootCI,label=label)))
} # end of bootmature

#' @title comparetwo compares two maturity curves fitted to two datasets
#' 
#' @description comparetwo takes two datasets selected on a categorical factor
#'    then fits a maturity curve mature = length + comparison, and compares 
#'    that with a maturity curve mature = length. This requires the MASS library
#'    for the version of anaova that permits test='Chisq', which conducts an
#'    analysis of deviance (see Venables and Ripley, 2nd ed, p288). The function
#'    returns the expanded model, the combined model, and a maturity model for 
#'    each of the separate data sets. 
#'
#' @param samdat1 the selected maturity data for the first set
#' @param samdat2 the selected maturity data for the second set
#' @param mature the name of the data.frame column containing the 1's and 0's
#'    denoting mature and immature, default = 'mature'
#' @param length the data.frame column name of the length variable, defualt=
#'    'length'
#' @param comparison the data.frame name of the comparison factor, default=
#'    'season' the only current alternative is 'site'
#' @param lower the smallest size class for inclusion in the predicted plot 
#'     data, the default = 50
#' @param upper the largest size class for inclusion in the predicted plot 
#'     data, the default = 160
#'
#' @return a list of lists containing all models, the anova, and the probabilty  
#'    the curves are different.
#' @export
#'
#' @examples
#' \dontrun{
#'   dat1 <- selecttwo(samdat=sam,value1=site,value2=c(seas1),
#'                     field1="site",field2="season")
#'   dat2 <- selecttwo(samdat=sam,value1=site,value2=c(seas2),#c("Autumn","Winter"),
#'                     field1="site",field2="season")
#'   ans <- comparetwo(dat1,dat2)
#'   str(ans,max.level=1)
#' }
comparetwo <- function(samdat1,samdat2,mature="mature",length="length",
                       comparison="season",lower=50,upper=160) {
  samdat <- rbind(samdat1,samdat2)
  samdat[,comparison] <- factor(samdat[,comparison])
  switch(comparison,
         "season" = modelS <-glm(mature ~ length + season,  data=samdat, 
                                 family=binomial(link="logit")),
         "site" = modelS <-glm(mature ~ length + site,  data=samdat, 
                               family=binomial(link="logit")),
         stop("Comparison not included"))
  propm <- tapply(samdat[,mature],samdat[,length],mean) #mean maturity at L
  lens <- as.numeric(names(propm)) # lengths in the data
  out <- binglm(modelS,console=FALSE)
  param <- out$param
  predm <- maturity(param[1],param[2],lens)
  plotdat <- cbind(lens,propm,predm)
  modelC <- list(out=out,model=modelS,param=param,plotdat=plotdat)
  
  model <- glm(mature ~ length, data=samdat, family = binomial(link="logit"))
  anv <- anova(modelS, model, test = "Chisq")
  prob <- anv[,"Pr(>Chi)"][2]
  model <- fitmaturity(samdat,mature=mature,length=length,lower=lower,upper=upper)
  first <- fitmaturity(samdat1,mature=mature,length=length,lower=lower,upper=upper)
  second <- fitmaturity(samdat2,mature=mature,length=length,lower=lower,upper=upper)
  return(invisible(list(model=model,modelC=modelC,anv=anv,prob=prob,
                        first=first,second=second))) 
} # end of comparetwo


#' @title comparetwogroups compares two groups of sites using logistic regression
#' 
#' @description comparetwogroups using logistic regressions to compare whether 
#'     two groups of sites are sufficiently similar for their data to be 
#'     combined, or not.
#'
#' @param grp1 a vector of at least one site number to be found in the 
#'     data.frame samdat. The default column name = 'site'.
#' @param grp2 a vector of at least one site number to be found in the 
#'     data.frame samdat. The default column name = 'site'.    
#' @param samdat a data.frame of maturity data with site numbers
#' @param sitecol column name of the variable for the site numbers. The default
#'     column name = 'site'
#' @param length column name in samdat for the length variable
#' @param mature column name for the maturity variable in samdat mature=1,
#'     immature = 0
#' @param sigdiff what probability level constistutes a significant difference.
#'     default=0.05
#' 
#' @seealso {
#'    \link{getloc}, \link{siteprops}, \link{binglm}
#' }
#'
#' @return a list of the anova between the models with and without group/site,
#'     the significance of any difference, model1, model2, and the outputs from
#'     applying binglm to each model
#' @export
#'
#' @examples
#' args(comparetwogroups)   
#' print("wait on internal data sets")
#' # grp1=sitenames[1];grp2=sitenames[2]; sitecol="siteseas";samdat=sam;
#' # length="length";mature="mature"; sigdiff=0.05
comparetwogroups <- function(grp1,grp2,samdat,sitecol="site",
                             length="length",mature="mature",sigdiff=0.05) { 
  picks <- which(samdat[,sitecol] %in% c(grp1,grp2))
  sams <- samdat[picks,]
  samsf <- sams
  samsf$sitenum <- NA
  picks1 <- which(samsf[,sitecol] %in% grp1)
  samsf$sitenum[picks1] <- 1
  samsf$sitenum[-picks1] <- 2
  samsf$sitenum <- factor(samsf$sitenum)
  model1 <- glm(formula = samsf[,mature] ~ samsf[,length] + samsf$sitenum, 
                family = binomial)
  model2 <- glm(formula = samsf[,mature] ~ samsf[,length],family = binomial)
  anv <- anova(model1, model2, test = "Chisq")
  prob <- anv[,"Pr(>Chi)"][2]
  signif <- FALSE
  if (prob < sigdiff) signif <- TRUE
  summ1 <- binglm(model1,console=FALSE)
  summ2 <- binglm(model2,console=FALSE)
  return(list(anv=anv,signif=signif,prob=prob,model1=model1,model2=model2,
              summ1=summ1,summ2=summ2))
} # end of comparetwogroups

#' @title findgroups sequentially compares sites and groups by maturity curve
#' 
#' @description findgroups conducts a sequential comparison of logistic maturity
#'     model fits and combines site numbers with statistically similar curves.
#'     Ideally, to test for homogeneity of the groups found, whatever ordering 
#'     is used shoulkd be reversed and the findgroups applied again.
#'
#' @param sitenums a vector of multiple site number from the data.frame samdat. 
#'     They should be ordered according to Lm50, latitude or longitude, or some 
#'     other factor that might suggest grouping
#' @param samdat a data.frame of maturity data with site numbers
#' @param sitecol column name of the variable for the site numbers. The default
#'     column name = 'site' 
#' @param console should the group structure be echoed to the console as they 
#'     are found? default=FALSE
#'
#' @return a list of groups of site numbers found to be statistically close 
#' @export
#'
#' @examples
#' args(findgroups)
#' print("wait on internal data sets")
#' # sitenums=sitenums;samdat=samb; sitecol="site";console=FALSE
findgroups <- function(sitenums,samdat,sitecol="site",console=FALSE) {
  nsites <- length(sitenums)
  groups <- NULL
  gn <- 1
  count <- 2
  grp <- sitenums[1]
  while (count <= nsites) {
    signif <- FALSE
    repeat {
      newsite <- sitenums[count]
      ans <- comparetwogroups(grp1=c(grp),grp2=newsite,samdat=samdat,
                             sitecol=sitecol)
      if (ans$signif) break
      grp <- c(grp,newsite)
      count <- count + 1
      if (count > nsites) break
    }
    groups[[gn]] <- grp
    if (console) print(groups[[gn]])
    if (count <= nsites) grp <- sitenums[count]
    gn <- gn + 1
    if (count == nsites) groups[[gn]] <- grp
    count <- count + 1
  }
  return(groups)
} # end of findgroups

#' @title fitgroups fits a maturity ~ length logistic regression to groups
#' 
#' @description fitgroups after using findgroups one can use fitgroups to 
#'     fit maturity curves to each group found 
#'
#' @param groups a list containing vectors of site numbers
#' @param samdat a data.frame of maturity data containing the required site(s)
#' @param sitecol column name of the variable for the site numbers. The default
#'     column name = 'site'
#' @param length column name in samdat for the length variable
#' @param mature column name for the maturity variable in samdat mature=1,
#'     immature = 0
#'
#' @return an invisible list containing the parameters and the fitmaturity 
#'     results
#' @export
#'
#' @examples
#' args(fitgroups)
#' print("wait on internal data sets")
#' #  groups=group; samdat=samb;length="length";mature="mature";sitecol="site"
fitgroups <- function(groups,samdat,sitecol="site",length="length",
                      mature="mature") {
  ng <- length(groups)
  fits <- vector(mode="list",length=ng)
  names(fits) <- 1:ng
  columns <- c("a","b","L50","IQ","n","nsites")
  param <- matrix(0,nrow=ng,ncol=length(columns),dimnames=list(1:ng,columns))
  for (i in 1:ng) { #  i=1
    pickg <- which(samdat[,sitecol] %in% groups[[i]])
    samg <- samdat[pickg,]
    ans <- fitmaturity(sams=samg,length=length,mature=mature)
    fits[[i]] <- ans
    param[i,] <- c(ans$param,(ans$out$out$df.null + 1),length(groups[[i]]))
  }
  return(invisible(list(param=param,fits=fits)))
} # end of fitgroups 

#' @title fitmaturity fits a maturity ~ length logistic regression to maturity
#' 
#' @description fitmaturity uses data of the form maturity = 1 immaturity = 0,
#'     and length and fits a glm using binomial errors (a logistic regression).
#'     It outputs the regression results, the parameters, a, b, Lm50, and IQ, as
#'     well as the data for plotting as 'lens', 'propm', and 'predm' 
#'
#' @param sams a data.frame of maturity data. This requires that any data 
#'     selection of say sites and/or seasons must be done outside the function. 
#' @param length column name in sams for the length variable
#' @param mature column name for the maturity variable in sams mature=1,
#'     immature = 0
#' @param lower the smallest size class for inclusion in the predicted plot 
#'     data, the default = 50
#' @param upper the largest size class for inclusion in the predicted plot 
#'     data, the default = 160
#'
#' @return an invisible list containing the regression results, the model,
#'     the parameters, and the plotdata, which contains both predm and 95 CI.
#' @export
#'
#' @examples
#' args(fitmaturity)
#' print("wait on internal data sets")
#' # sams=samdat;length="length";mature="mature";lower=50;upper=160
fitmaturity <- function(sams,length="length",mature="mature",lower=50,upper=160) {
  model <- glm(mature ~ length, data=sams, family = binomial(link="logit"))
  propm <- tapply(sams[,mature],sams[,length],mean) #mean maturity at L
  lens <- as.numeric(names(propm)) # lengths in the data
  out <- binglm(model,console=FALSE)
  param <- out$param
  newdat <- data.frame(length=(lower:upper))
  predm <- predict(model, newdata=newdat, se.fit=TRUE)
  nplot <- length(predm$fit)
  propm2 <- rep(NA,nplot)
  picklens <- which((lens <= upper) & (lens >= lower))
  pickL <- match(lens[picklens],newdat[,"length"])
  propm2[pickL] <- propm[picklens]
  predmat <- exp(predm$fit)/(1+exp(predm$fit))
  L95 <- exp(predm$fit-1.96*predm$se.fit)/(1+exp(predm$fit-1.96*predm$se.fit))
  U95 <- exp(predm$fit+1.96*predm$se.fit)/(1+exp(predm$fit+1.96*predm$se.fit))
  plotdat <- cbind(newdat=newdat,propm=propm2,predm=predmat,se=predm$se.fit,
                   L95=L95,U95=U95)
  return(invisible(list(out=out,model=model,param=param,plotdat=plotdat)))
} # end of fitmaturity


#' @title getdistances estimates the haversine distance between sites
#' 
#' @description getdistance uses theoutput from siteprops to estimate the
#'     haversine distances between sites. It output a symmetrical matrix of
#'     these values.
#'
#' @param siteprop the output from 'siteprops' containing 'combsite', 'lat' and
#'     'long' as a minimum
#'
#' @return a square symmetric matrix of haversine distances (kms) between all
#'     sites in a block
#' @export
#'
#' @examples
#' args(getdistances)
#' print("wait on internal data sets")
getdistances <- function(siteprop) {
  nsites <- nrow(siteprop)
  dist <- matrix(0,nrow=nsites,ncol=nsites,
                 dimnames=list(siteprop$site,siteprop$site))
  for (i in 1:nsites) {  # i = 1; j=2
    for (j in 1:nsites) {
      dist[j,i] <- as.numeric(haversine(siteprop[i,"lat"],siteprop[i,"long"],
                                        siteprop[j,"lat"],siteprop[j,"long"]))[2]
    }
  }
  return(dist)
} # end of getdistances

#' @title getloc pulls out the lat longs and block for identified site numbers
#' 
#' @description getloc isa utility function that extracts the lat longs, if 
#'     present, for all the identified site numbers
#'
#' @param sits a vector of site numbers to be found in the data.frame samdat
#' @param samdat a data.frame of maturity data with site numbers
#' @param orderby this can be 'lat', 'long', or 'none'. If 'none' no reordering
#'     will be done, otherwise the sites will be sorted by lat or long. default
#'     = 'none'
#'
#' @return a matrix of site by lat, long, block, and sitename
#' @export
#'
#' @examples
#' args(getloc)
#' print("wait on internal data sets")
#' # sits=bsites; samdat=samb; orderby="none"
getloc <- function(sits,samdat,orderby="none") { # 
  nsites <- length(sits)
  columns <- c("lat","long","block","sitename","combsite")
  numcol <- length(columns)
  locations <- as.data.frame(matrix(0,nrow=nsites,ncol=numcol,
                                    dimnames=list(sits,columns)))
  for (i in 1:nsites) {
    pick <- which(samdat$site == sits[i])
    locations[i,] <- c(samdat[pick[1],"lat"],samdat[pick[1],"long"],
                       samdat[pick[1],"block"],samdat[pick[1],"sitename"],
                       samdat[pick[1],"site"]) 
  }
  if (orderby %in% c("lat","long"))
    locations <- locations[order(locations[,orderby]),]
  return(locations)
} # end of getloc

#' @title getsiteparams get all the individual site parameters in a data.frame
#' 
#' @description getsiteparams estimates the site parameters and all the model
#'     fits for all sites within a data.frame. This sets up the system ready
#'     for site comparisons 
#'
#' @param samdat a data.frame of maturity data containing sites
#' @param sitecol column name of the variable for the site numbers. The default
#'     column name = 'site'
#' @param length column name in samdat for the length variable
#' @param mature column name for the maturity variable in samdat mature=1,
#'     immature = 0
#' @param sau column name for the sau or block in samdat, default='block'
#' @param lower the smallest size class for inclusion in the predicted plot 
#'     data, the default = 50
#' @param upper the largest size class for inclusion in the predicted plot 
#'     data, the default = 160
#'
#' @return a list containing a matrix of the site parameters and a list of all
#'     the model fits used (called allfit)
#' @export
#'
#' @examples
#' args(getsiteparams)
#' print("wait on internal data sets")
#' samdat=samnw; sitecol="site";length="length";mature="mature";sau="block";lower=20;upper=160
getsiteparams <- function(samdat,sitecol="site",length="length",mature="mature",
                          sau="block",lower=50,upper=160) {
  sits <- sort(unique(samdat[,sitecol]))
  nsits <- length(sits)
  allfit <- vector(mode="list",length=nsits)
  names(allfit) <- sits
  columns <- c("a","b","Lm50","IQm","n","site","sau")
  sitepar <- matrix(0,nrow=nsits,ncol=length(columns),dimnames=list(sits,columns))
  for (i in 1:nsits) { # i = 1
    picks <- which(samdat[,sitecol] == sits[i])
    model <- fitmaturity(sams=samdat[picks,],length=length,mature=mature,
                         lower=lower,upper=upper)
    allfit[[i]] <- model
    sitepar[i,] <- c(model$param,length(model$out$out$deviance.resid),
                                        sits[i],NA)
  }
  area <- sort(unique(samdat[,sau]))
  narea <- length(area)  
  if (narea > 0) {
    if (narea == 1) {
      sitepar[,"sau"] <- area
    } else {
     sitebysau <- table(samdat[,sitecol],samdat[,sau])
     sits <- as.numeric(rownames(sitebysau))
     for (i in 1:narea) { # i = 1
       pick <- which(sitebysau[,as.character(area[i])] > 0)
       sitepar[pick,"sau"] <- area[i]
     }
   }
  }
  return(invisible(list(sitepar=sitepar,allfit=allfit)))
} # end of getsiteparams

#' @title groupprops tabulates the properties of each site number
#' 
#' @description groupprops extracts the main properties of two vectors of site
#'     numbers. These properties include: n, site, sitename, lat, long, block,
#'     subblock, date, and propM (proportion mature). It uses siteprops
#'
#' @param grp a vector of at least one site number to be found in the 
#'     data.frame samdat. The default column name = 'site'.
#' @param samdat a data.frame of maturity data with site numbers
#' @param groupname a name for the group, default='grp'
#' @param sitecol column name of the variable for the site numbers. The default
#'     column name = 'site'
#' 
#' @seealso {
#'    \link{getloc}, \link{findgroups}, \link{fitgroups}
#' }
#'
#' @return a list of the sites, the block(s), n, propI and propM
#' @export
#'
#' @examples
#' args(groupprops)
#' print("wait on internal data sets")
#' # grp=c(780,783); samdat=samw; sitecol="site"; groupname="blk13"
groupprops <- function(grp,samdat,groupname="grp",sitecol="site") {
  picks <- which(samdat[,sitecol] %in% grp)
  tot <- length(picks)
  pickM <- which(samdat$mature[picks] == 1)
  if (length(pickM) == 0) {
    propI <- 1
    propM <- 0
  } else {
    propI <- (tot - length(pickM))/tot
    propM <- length(pickM)/tot
  }
  blks <- sort(unique(samdat$block[picks]))
  
  if (length(blks) > 1) blks <- makelabel("",blks)
  columns <- c("sites","block","n","propI","propM")  
  ans <- as.data.frame(matrix(0,nrow=1,ncol=length(columns),
                              dimnames=list(groupname,columns)))
  ans[1,] <- c(makelabel("",grp),blks,tot,propI,propM)
  return(ans)
} # end of groupprops

#' @title labelmake combines the contents of a character vector to make a label
#' 
#' @description labelmake is a utility function that takes a vector of values, 
#'     either numeric or text, and combines them with a selected separator. Any
#'     numeric values also have the option of being rounded to a given precision
#'
#' @param vect a vector of either all text or all numbers, mixtures not allowed
#' @param sep what separator should be used, default = '_'
#' @param digits if a numeric vector is input how many decimal places to be used?
#'     default = 3
#'
#' @return a character string of the contents of vect combined
#' @export
#'
#' @examples
#' x <- c("one","two","three")
#' labelmake(x)
#' x <- c(1,2,3,pi)
#' labelmake(x,digits=5)
labelmake <- function(vect,sep="_",digits=3) {
  if (is.numeric(vect)) vect <- round(vect,digits)
  label <- paste(vect,sep=sep,collapse=sep)
  return(label)
} # end of labelmake

#' @title maturity Logistic maturity curve
#'
#' @description maturity this uses the logistic function:
#'   exp(a+bL)/(1+exp(a+bL)), which has the property that the SM50 = -a/b
#'   and the interquartile distance is 2.Ln(3)/b.
#' @param ina is the intercept of the exponential function
#' @param inb is the gradient of the exponential function
#' @param lens a vector of lengths for which the logistic maturity value
#'     will be calculated
#' @return A vector of length(lens) containing the predicted maturity at
#'     length values
#' @export
#'
#' @examples
#' a <- -14.383
#' b <- 0.146017
#' lens <- seq(2,210,2)
#' Maturity <- maturity(a,b,lens)
maturity <- function(ina,inb,lens) {
  ans <- exp(ina+inb*lens)/(1+exp(ina+inb*lens))
  return(ans)
} # end of maturity

#' @title plotcomparison generates a plot from the output of comparetwo
#' 
#' @description plotcomparison uses the output of the function 'comparetwo' to
#'    generate a plot of the maturity curves fitted to the original data from 
#'    the two datasets, plus, if there is no significant difference, the 
#'    combined maturity curve, all with 95 percent CI.
#'
#' @param ans the output from the function 'comparetwo'
#' @param name1 a character strong as the name of the first dataset
#' @param name2 a character strong as the name of the second dataset
#' @param P the probabiolity lelve that is taken as representing a significant
#'    difference between two maturity curves, default = 0.01
#'
#' @return nothing but it does plot a graph
#' @export
#'
#' @examples
#' print("wait on internal data sets")
plotcomparison <- function(ans,name1, name2,P=0.01) {
  first <- ans$first
  plotdat1 <- first$plotdat
  second <- ans$second
  plotdat2 <- second$plotdat
  model <- ans$model
  plotdat <- model$plotdat
  lens <- plotdat2[,"length"]
  parset(cex=1.0)
  plot(lens,plotdat1[,"propm"], type="p",pch=16,cex=1.0, ylim=c(0, 1),col=1, 
       ylab="Proportion Mature", xlab="length",panel.first=grid(),
       xlim=c(50,160))
  lines(lens, plotdat1$predm, col="black",lwd=2)
  lines(lens, plotdat1$L95, lty=2)
  lines(lens, plotdat1$U95, lty=2)
  
  points(lens,plotdat2[,"propm"],pch=1,cex=1.5,col=4)
  lines(lens, plotdat2$predm, col=4,lwd=2)
  lines(lens, plotdat2$L95, lty=2,col=4)
  lines(lens, plotdat2$U95, lty=2,col=4)
  label <- c(name1,name2)
  linecol <- c(1,4)
  if (ans$prob > P) {
    lines(lens, plotdat$predm, col=2,lwd=2)
    lines(lens, plotdat$L95, lty=2,col=2)
    lines(lens, plotdat$U95, lty=2,col=2)
    label <- c(name1,name2,"combined")
    linecol <- c(1,4,2)
  }
  abline(h=0.5,lwd=1,col="darkgrey")
  legend("topleft",legend=label,col=linecol,lwd=3,bty="n")
  plab <- paste0("pdiff = ",round(ans$prob,4),"    ")
  mtext(plab,side=1,outer=FALSE,line=-1.1,adj=1)
} # end of plotcomparison


#' @title plotgroups fits a maturity ~ length logistic regression to groups
#' 
#' @description plotgroups after using findgroups one can use plotgroups to 
#'     plot maturity curves to each group found. It uses fitgroups and outputs
#'     those results invisibly.
#'
#' @param groups a list containing vectors of site numbers
#' @param samdat a data.frame of maturity data containing the required site(s)
#' @param xmin a generic x-axis minimum?, if 0 use local, default=0
#' @param xmax a generic x-axis maximum?, if 0 use local, default=0
#' @param sitecol column name of the variable for the site numbers. The default
#'     column name = 'site'
#' @param length column name in samdat for the length variable
#' @param mature column name for the maturity variable in samdat mature=1,
#'     immature = 0
#' @param margin what margin should be used in the plot? 
#'     default=c(0.4,0.4,0.05,0.05)
#'
#' @return an invisible list containing the parameters and the fitmaturity 
#'     results. It also generates a plot.
#' @export
#'
#' @examples
#' args(plotgroups)
#' print("wait on internal data sets")
#' #  groups=group; samdat=samb;length="length";mature="mature";sitecol="site"
plotgroups <- function(groups,samdat,xmin=0,xmax=0,sitecol="site",
                       length="length",mature="mature",
                       margin=c(0.4,0.4,0.05,0.05)) { # x = outgs
  outfit <- suppressWarnings(fitgroups(groups,samdat,sitecol=sitecol,
                                       length=length,mature=mature))
  fits <- outfit$fits
  nfits <- length(fits)
  parset(plots=pickbound(nfits),byrow=FALSE,margin=margin)
  for (i in 1:nfits) {
    ylabel <- makelabel("",groups[[i]])
    ans <- plotmaturity(fits[[i]],label=ylabel,xmin=xmin,xmax=xmax)
  }
  return(invisible(outfit))
} # end of plotgroups

#' @title plotmaturity generates a plot of fitted maturity data.
#' 
#' @description plotmaturity generates a plot of the maturity curve and data 
#'    when provided with the model output from fitmaturity
#'
#' @param model the output from the fitmaturity function
#' @param label the label for the y-axis, usually the site number or the
#'     combination of site numbers making up a group
#' @param col colour of the predicted maturity curve, default=2 = red
#' @param xmin a generic x-axis minimum?, if 0 use local, default=0
#' @param xmax a generic x-axis maximum?, if 0 use local, default=0
#' @param CI should the 95 perc CI be plotted? default=FALSE
#' @param setpar should parset be called, useful when only a single plot is
#'     being made. default=FALSE
#' @param lwd what line width to use for predicted maturity line, default=2
#' 
#' @seealso {
#'   \link{fitmaturity}, \link{plotsinglegroup} 
#' }
#'
#' @return nothing but it does generate a plot
#' @export
#'
#' @examples
#' args(plotmaturity)
#' print("wait on data sets")
plotmaturity <- function(model,label="",col=2,xmin=0,xmax=0,CI=FALSE,
                         setpar=FALSE,lwd=2) {
  dat <- model$plotdat
  nobs <- length(model$model$residuals)
  xlabel <- ""
  if (setpar) {
    parset()
    xlabel <- "Shell Length mm"
  }
  if ((xmin == 0) & (xmax == 0)) {
    xmin <- getmin(dat[,"length"])
    xmax <- getmax(dat[,"length"])
  } 
  plot(dat[,"length"],dat[,"propm"],type="p",cex=1.0,pch=16,xlim=c(xmin,xmax),
       xlab=xlabel,ylab=label,yaxs="r",panel.first = grid())
  lines(dat[,"length"],dat[,"predm"],lwd=lwd,col=col)
  abline(h=c(0,0.5,1),lwd=1,col="grey")
  Lm50 <- model$param[1,"L50"]
  abline(v=Lm50,lwd=1,col=3)
  if (CI) {
    lines(dat[,"length"],dat[,"U95"],lwd=1,col=1,lty=2)
    lines(dat[,"length"],dat[,"L95"],lwd=1,col=1,lty=2)
  }
  mtext(round(Lm50,3),side=3,line=-1.2,cex=1.0,adj=0)
  mtext(nobs,side=3,line=-2.5,cex=1.0,adj=0)
} # end of plotmaturity

#' @title plotsinglegroup plots a maturity ~ length logistic regression to groups
#' 
#' @description plotsinglegroup after using findgroups one can use plotgroup to 
#'     plot maturity curves to each site in a group alongside an overal model 
#'     fit.
#'
#' @param groups a vector of site numbers
#' @param samdat a data.frame of maturity data containing the required site(s)
#' @param xmin a generic x-axis minimum?, if 0 use local, default=0
#' @param xmax a generic x-axis maximum?, if 0 use local, default=0
#' @param sitecol column name of the variable for the site numbers. The default
#'     column name = 'site'
#' @param length column name in samdat for the length variable
#' @param mature column name for the maturity variable in samdat mature=1,
#'     immature = 0
#' @param lwd what line width to use for predicted maturity line, default=2   
#' @param makeplot should a plot be produced, default=TRUE  
#'
#' @return an invisible list containing the parameters and the fitmaturity 
#'     results. It also generates a plot.
#' @export
#'
#' @examples
#' args(plotsinglegroup)
#' print("wait on internal data sets")
#' #  groups=sits; samdat=sams;length="length";mature="mature";sitecol="site";xmin=0;xmax=0
plotsinglegroup <- function(groups,samdat,xmin=0,xmax=0,sitecol="site",
                            length="length",mature="mature",lwd=2,makeplot=TRUE) { 
  outfit <- fitgroups(groups,samdat,sitecol=sitecol,length=length,mature=mature)
  fits <- outfit$fits
  nfits <- length(fits)
  combined <- fitmaturity(sams=samdat,length=length,mature=mature)
  if (makeplot) {
    label <- paste0("sau",samdat$block[1])
    plotmaturity(combined,label=label,xmin=xmin,xmax=xmax,lwd=lwd)
    colour <- 3:(nfits+2)
    for (i in 1:nfits) { #  i = 1
      plotdat <- fits[[i]]$plotdat
      lines(plotdat[,"length"],plotdat[,"predm"],lwd=1,col=colour[i])
    }
    legend("bottomright",legend=c(groups,nfits),col=c(colour,2),lwd=3,
           bty="n",cex=1.25)
  }
  return(invisible(list(outfit=outfit,combined=combined)))
} # end of plotsinglegroup

#' @title plottwosites fits and plots two sites, including the combined by default
#' 
#' @description plottwosites fits a maturity curve for two sites and, by default, 
#'     the combined data for both sites. It then plots them together on the same 
#'     graph. The out is, invisibly, a list of all results.
#'
#' @param sam1 a data.frame containing the maturity data for group 1
#' @param sam2 a data.frame containing the maturity data for group 2
#' @param length column name of the length variable within the data.frame,
#'     default='length'
#' @param mature column name of the maturity variable containing mature=1, and
#'     immature = 0 for all length observations. default='mature'
#' @param sitecol column name of the variable for the site numbers. The default
#'     column name = 'site'
#'
#' @return invisibly, a list containing the model fits for each group of sites, 
#'     and the model for all sites combined
#' @export
#'
#' @examples
#' args(plottwogroups)
#' print("wait on internal data sets")
#' # grp1=c(43,874,2);grp2=c(128); samdat=samW;length="length";mature="mature";sitecol="site"
plottwogroups <- function(sam1,sam2,length="length",mature="mature",
                          sitecol="site") {  
  sams <- rbind(sam1,sam2)
  ans <- vector("list",length=3)
  ans[[1]] <- fitmaturity(sam1,length=length,mature=mature)
  ans[[2]] <- fitmaturity(sam2,length=length,mature=mature)
  ans[[3]] <- fitmaturity(sams,length=length,mature=mature)
  xmin <- getmin(c(ans[[1]]$plotdat[,"length"],ans[[2]]$plotdat[,"length"]))
  xmax <- getmax(c(ans[[1]]$plotdat[,"length"],ans[[2]]$plotdat[,"length"]))
  parset()
  plot(ans[[1]]$plotdat[,"length"],ans[[1]]$plotdat[,"propm"],type="p",pch=16,
       cex=1.0,xlim=c(xmin,xmax),xlab="shell length",yaxs="r",
       ylab=paste0("Proportion Mature"),panel.first=grid())
  points(ans[[2]]$plotdat[,"length"],ans[[2]]$plotdat[,"propm"],pch=1,
         col="blue",cex=1.75)
  lines(ans[[1]]$plotdat[,"length"],ans[[1]]$plotdat[,"predm"],lwd=2,col=1)
  lines(ans[[2]]$plotdat[,"length"],ans[[2]]$plotdat[,"predm"],lwd=2,col="blue")
  points(ans[[3]]$plotdat[,"length"],ans[[3]]$plotdat[,"propm"],pch=1,col=2,
         cex=1.5)
  lines(ans[[3]]$plotdat[,"length"],ans[[3]]$plotdat[,"predm"],lwd=2,col=2)
  lab1 <- makelabel("",sort(unique(sam1[,sitecol])))
  lab2 <- makelabel("",sort(unique(sam2[,sitecol])))
  legend("topleft",c(lab1,lab2,"combined"),col=c("black","blue",2),
         lwd=3,bty="n",cex=1.2)
  param <- sapply(ans,"[[","param")
  rownames(param) <- colnames(ans[[1]]$param)
  colnames(param) <- c(lab1,lab2,"combined")
  ans$param <- param
  return(invisible(ans))
} # end of plottwogroups

#' @title printgroup a simple function to print the group list 
#' 
#' @description printgroup is a simple function to print a group of sites list
#'     that is produced by the function 'findgroups'. This is used instead of
#'     defining a new S3 method so one doesn't have to define a new class.  
#'
#' @param group a list of groups of sites from the maturity data. Usually
#'     derived from the function 'findgroups'
#'     
#' @seealso{
#'   \link{findgroups}
#' }
#'
#' @return nothing but it does print out the group list more tidily
#' @export
#'
#' @examples
#' args(printgroup)
#' agroup <- list("1"=c(1,2,3),"2"=c(4,5),"3"=c(6,7,8,9))
#' printgroup(agroup)
printgroup <- function(group) {
  if (is.null(group)) {
    print("NULL")
  } else {
    ngrp <- length(group)
    for (i in 1:ngrp) {
      cat(paste0("group_",i),":  ",sort(group[[i]]),"\n")
    }
  }
} # end of printgroup

#' @title selecttwo selects data from a data.frame on two factors
#' 
#' @description selecttwo selects data from a data.frame on two factors/fields
#'     where the values for each field can be a vector of values. This can be 
#'     useful when there ar emultiple levels within multiple fields/factors and
#'     how best to combine them is unknown.
#'
#' @param samdat the data.frame being selected from.
#' @param value1 a vector or values for the first field, possibly just one
#' @param value2 a vector or values for the second field, possibly just one
#' @param field1 the column name in the data.frame from which to select the
#'     value1
#' @param field2 the column name in the data.frame from which to select the
#'     value2
#'
#' @return a data.frame subsetted from the input data.rame
#' @export
#'
#' @examples
#' \dontrun{
#'   site1 <- "George"
#'   seas1 <- c("Spring","Summer")
#'   dat1 <- selecttwo(samdat=sam,value1=site1,value2=c(seas1),
#'                     field1="site",field2="season")
#' }
selecttwo <- function(samdat,value1,value2,field1="site",field2="season") {
  pickrow <- which((samdat[,field1] %in% value1) & (samdat[,field2] %in% value2))
  outdat <- samdat[pickrow,]
  return(outdat)                 
} # end of selecttwo

#' @title siteprops tabulates the properties of each site number in a data.frame
#' 
#' @description siteprops extracts the main properties of a vector of site
#'     numbers. These properties include: n, site, sitename, lat, long, block,
#'     subblock, date, and propM (proportion mature).
#'
#' @param samdat a data.frame of maturity data with site numbers
#' @param sitecol column name of the variable for the site numbers. The default
#'     column name = 'site'
#' 
#' @seealso {
#'    \link{getloc}, \link{comparetwogroups}
#' } 
#'
#' @return a data.frame of site properties
#' @export
#'
#' @examples
#' args(siteprops)
#' print("wait on internal data sets")
#' # sits=bsites; samdat=samb
siteprops <- function(samdat,sitecol="site") {
  sits <- sort(unique(samdat[,sitecol]))
  nsites <- length(sits)
  columns <- c("n","site","sitename","lat","long","block","subblock","date",
               "propM")
  ans <- as.data.frame(matrix(0,nrow=nsites,ncol=length(columns),
                              dimnames=list(sits,columns)))
  for (i in 1:nsites) {  # i = 1
    picks <- which(samdat[,sitecol] == sits[i])
    sams <- samdat[picks,]
    n <- length(picks)
    ans[i,] <- c(n,sits[i],sams[1,"sitename"],sams[1,"lat"],sams[1,"long"],
                 sams[1,"block"],sams[1,"subblock"],sams[1,"date"],sams[1,"propM"])
  }
  label <- c("n","lat","long","block","propM")
  for (i in 1:5) ans[,label[i]] <- as.numeric(ans[,label[i]])
  return(ans)
} # end of siteprops

#' @title testgroup literally compares sites using logistic regression
#' 
#' @description testgroup using logistic regressions to compare whether 
#'     two or more sites are sufficiently similar for their data to be combined, 
#'     or not. 
#'
#' @param grp a vector of at least two site numbers to be found in the 
#'     data.frame samdat. The default column name = 'site'.
#' @param samdat a data.frame of maturity data with site numbers
#' @param sitecol column name of the variable for the site numbers. The default
#'     column name = 'site'
#' @param length column name in samdat for the length variable
#' @param mature column name for the maturity variable in samdat mature=1,
#'     immature = 0
#' 
#' @seealso {
#'    \link{getloc}, \link{siteprops}
#' }
#'
#' @return a list of a summary of the logistic regression, and the model
#' @export
#'
#' @examples
#' args(testgroup)
#' print("wait on internal data sets")
#' # site1=grp;site2=grp2; sitecol="site";samdat=samb;length="length";mature="mature"
testgroup <- function(grp,samdat,sitecol="site",length="length",mature="mature") { 
  picks <- which(samdat[,sitecol] %in% grp)
  sams <- samdat[picks,]
  samsf <- sams
  samsf[,sitecol] <- factor(samsf[,sitecol])
  model <- glm(formula = samsf[,mature] ~ samsf[,length] + samsf[,sitecol], 
               family = binomial)
  summ <- summary(model)
  return(summ=summ)
} # end of testgroup
