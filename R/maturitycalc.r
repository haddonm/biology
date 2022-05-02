
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
#' print("wait on internal datasets")
#' # blk=13; samdat=samw; sitecol="site"
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
#' \dontrun{
#'  # wait on a data set being included into the package
#' }
binglm <- function(x,digits=6,console=TRUE) { # x=smodel; digits=6;console=TRUE
  out <- summary(x)
  pars <- coef(out)[,"Estimate"]
  npar <- length(pars)
  nmod <- npar-1
  columns <- c("a","b","L50","IQ")
  rows <- names(pars[1:nmod])
  models <- matrix(0,nrow=nmod,ncol=4,dimnames=list(rows,columns))
  models[1,] <- c(pars[1],pars[npar],-pars[1]/pars[npar],2*log(3)/pars[npar])
  if (nmod > 1) {
    for (i in 2:nmod) {
      par1 <- pars[1] + pars[i]
      models[i,] <- c(par1,pars[npar],-par1/pars[npar],2*log(3)/pars[npar]) 
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
#' @param mat10 column name for the maturity variable in samdat mature=1,
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
#' print("wait on internal data sets")
#' # grp1=g1;grp2=g2; sitecol="site";samdat=samb;
#' # length="length";mat10="mature"; sigdiff=0.05
comparetwogroups <- function(grp1,grp2,samdat,sitecol="site",
                             length="length",mat10="mature",sigdiff=0.05) { 
  picks <- which(samdat[,sitecol] %in% c(grp1,grp2))
  sams <- samdat[picks,]
  samsf <- sams
  samsf$sitenum <- NA
  picks1 <- which(samsf$site %in% grp1)
  samsf$sitenum[picks1] <- 1
  samsf$sitenum[-picks1] <- 2
  samsf$sitenum <- factor(samsf$sitenum)
  model1 <- glm(formula = samsf[,mat10] ~ samsf[,length] + samsf$sitenum, 
                family = binomial)
  model2 <- glm(formula = samsf[,mat10] ~ samsf[,length],family = binomial)
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
#' @param mat10 column name for the maturity variable in samdat mature=1,
#'     immature = 0
#'
#' @return an invisible list containing the parameters and the fitmaturity 
#'     results
#' @export
#'
#' @examples
#' print("wait on internal data sets")
#' #  groups=groups1; samdat=samb;length="length";mat10="mature";sitecol="site"
fitgroups <- function(groups,samdat,sitecol="site",length="length",
                      mat10="mature") {
  ng <- length(groups)
  fits <- vector(mode="list",length=ng)
  names(fits) <- 1:ng
  columns <- c("a","b","L50","IQ","n","nsites")
  param <- matrix(0,nrow=ng,ncol=length(columns),dimnames=list(1:ng,columns))
  for (i in 1:ng) { #  i=1
    ans <- fitmaturity(groups[[i]],samdat=samdat,sitecol=sitecol,
                       length=length,mat10=mat10)
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
#' @param sites a vector of at least one site number to be found in the 
#'     data.frame samdat. The default column name = 'site'.
#' @param samdat a data.frame of maturity data containing the required site(s)
#' @param sitecol column name of the variable for the site numbers. The default
#'     column name = 'site'
#' @param length column name in samdat for the length variable
#' @param mat10 column name for the maturity variable in samdat mature=1,
#'     immature = 0
#'
#' @return an invisible list containing the regression results, the model,
#'     the parameters, and the plotdata
#' @export
#'
#' @examples
#' print("wait on internal data sets")
#' #  sites=grp; samdat=samb;length="length";mat10="mature";sitecol="site"
fitmaturity <- function(sites,samdat,sitecol="site",length="length",
                        mat10="mature") {
  picks <- which(samdat[,sitecol] %in% sites)
  sams <- samdat[picks,]
  model <- glm(formula = sams[,mat10] ~ sams[,length], family = binomial)
  propm <- tapply(sams[,mat10],sams[,length],mean) #mean maturity at L
  lens <- as.numeric(names(propm)) # lengths in the data
  out <- binglm(model,console=FALSE)
  param <- out$param
  predm <- maturity(param[1],param[2],lens)
  plotdat <- cbind(lens,propm,predm)
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
#' @description getsiteparams estimates the site parmaeters and all the model
#'     fits for all sites within a data.frame. This sets up the system ready
#'     for site comparisons 
#'
#' @param samdat a data.frame of maturity data containing sites
#' @param sitecol column name of the variable for the site numbers. The default
#'     column name = 'site'
#' @param length column name in samdat for the length variable
#' @param mat10 column name for the maturity variable in samdat mature=1,
#'     immature = 0
#'
#' @return a list containing a matrix of the site parameters and a list of all
#'     the model fits used (called allfit)
#' @export
#'
#' @examples
#' print("wait on intenral data sets")
getsiteparams <- function(samdat,sitecol="site",length="length",mat10="mature") {
  sits <- sort(unique(samdat[,sitecol]))
  nsits <- length(sits)
  allfit <- vector(mode="list",length=nsits)
  names(allfit) <- sits
  columns <- c("a","b","Lm50","IQm","n")
  sitepar <- matrix(0,nrow=nsits,ncol=length(columns),dimnames=list(sits,columns))
  for (i in 1:nsits) { # i = 1
    model <- fitmaturity(sits[i],samdat=samdat,sitecol=sitecol,length=length,
                         mat10=mat10)
    allfit[[i]] <- model
    sitepar[i,] <- c(model$param,length(model$out$out$deviance.resid))
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
#' \dontrun{
#' a <- -14.383
#' b <- 0.146017
#' lens <- seq(2,210,2)
#' Maturity <- maturity(a,b,lens)
#' }
maturity <- function(ina,inb,lens) {
  ans <- exp(ina+inb*lens)/(1+exp(ina+inb*lens))
  return(ans)
} # end of maturity

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
#' @param mat10 column name for the maturity variable in samdat mature=1,
#'     immature = 0
#' @param margin what margin should be used in the plot? 
#'     default=c(0.4,0.4,0.05,0.05)
#'
#' @return an invisible list containing the parameters and the fitmaturity 
#'     results. It also generates a plot.
#' @export
#'
#' @examples
#' print("wait on internal data sets")
#' #  groups=group2; samdat=samb;length="length";mat10="mature";sitecol="site"
plotgroups <- function(groups,samdat,xmin=0,xmax=0,sitecol="site",
                       length="length",mat10="mature",
                       margin=c(0.4,0.4,0.05,0.05)) { # x = outgs
  outfit <- suppressWarnings(fitgroups(groups,samdat,sitecol=sitecol,
                                       length=length,mat10=mat10))
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
#' @param setpar should parset be called, useful when only a single plot is
#'     being made. default=FALSE
#' 
#' @seealso {
#'   \link{fitmaturity} 
#' }
#'
#' @return nothing but it does generate a plot
#' @export
#'
#' @examples
#' print("wait on data sets")
plotmaturity <- function(model,label="",col=2,xmin=0,xmax=0,setpar=FALSE) {
  dat <- model$plotdat
  xlabel <- ""
  if (setpar) {
    parset()
    xlabel <- "Shell Length mm"
  }
  if ((xmin == 0) & (xmax == 0)) {
    xmin <- getmin(dat[,"lens"])
    xmax <- getmax(dat[,"lens"])
  } 
  plot(dat[,"lens"],dat[,"propm"],type="p",cex=1.0,pch=16,xlim=c(xmin,xmax),
       xlab=xlabel,ylab=label,yaxs="r",panel.first = grid())
  lines(dat[,"lens"],dat[,"predm"],lwd=2,col=col)
  abline(h=c(0,0.5,1),lwd=1,col="grey")
  Lm50 <- model$param[1,"L50"]
  abline(v=Lm50,lwd=1,col=3)
  mtext(round(Lm50,3),side=3,line=-1.2,cex=1.0,adj=0)
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
#' @param margin what margin should be used in the plot? 
#'     default=c(0.4,0.4,0.05,0.05)
#' @param sitecol column name of the variable for the site numbers. The default
#'     column name = 'site'
#' @param length column name in samdat for the length variable
#' @param mat10 column name for the maturity variable in samdat mature=1,
#'     immature = 0
#'
#' @return an invisible list containing the parameters and the fitmaturity 
#'     results. It also generates a plot.
#' @export
#'
#' @examples
#' print("wait on internal data sets")
#' #  groups=groups1; samdat=samb;length="length";mat10="mature";sitecol="site"
plotsinglegroup <- function(groups,samdat,xmin=0,xmax=0,margin=c(0.4,0.4,0.05,0.05),
                            sitecol="site",length="length",mat10="mature") { 
  outfit <- fitgroups(groups,samdat,sitecol=sitecol,length=length,mat10=mat10)
  fits <- outfit$fits
  nfits <- length(fits)
  combined <- fitmaturity(sites=groups,samdat=samdat,sitecol=sitecol,
                          length=length,mat10=mat10) 
  label <- makelabel("",groups)
  parset(margin=margin)
  plotmaturity(combined,label=label,xmin=xmin,xmax=xmax)
  colour <- 3:(nfits+2)
  for (i in 1:nfits) { #  i = 1
    plotdat <- fits[[i]]$plotdat
    lines(plotdat[,"lens"],plotdat[,"predm"],lwd=2,col=colour[i])
  }
  legend("bottomright",c(groups,label),col=c(colour,2),lwd=3,bty="n",cex=1.25)
  return(invisible(list(outfit=outfit,combined=combined)))
} # end of plotsinglegroup

#' @title plottwosites fits and plots two sites, including the combined by default
#' 
#' @description plottwosites fits a maturity curve for two sites and, by default, 
#'     the combined data for both sites. It then plots them together on the same 
#'     graph. The out is, invisibly, a list of all results.
#'
#' @param grp1 a vector of at least one site number to be found in the 
#'     data.frame samdat. The default column name = 'site'.
#' @param grp2 a vector of at least one site number to be found in the 
#'     data.frame samdat. The default column name = 'site'.  
#' @param samdat a data.frame containing the maturity data that includes the two 
#'     sites
#' @param length column name of the length variable within the data.frame,
#'     default='length'
#' @param mat10 column name of the maturity variable containing mature=1, and
#'     immature = 0 for all length observations. default='mature'
#' @param sitecol column name for the site number, default='combsite', ready for 
#'     when sites are combined
#'
#' @return invisibly, a list containing the model fits for each site and the
#'     combined sites if combine = TRUE
#' @export
#'
#' @examples
#' print("wait on internal data sets")
#' # grp1=c(43,874,2);grp2=c(128); samdat=samW;length="length";mat10="mature";sitecol="site"
plottwogroups <- function(grp1,grp2,samdat,length="length",mat10="mature",
                         sitecol="site") {  
  picks <- which(samdat[,sitecol] %in% c(grp1,grp2))
  sams <- samdat[picks,]
  ans <- vector("list",length=3)
  ans[[1]] <- fitmaturity(grp1,sams,length=length,mat10=mat10)
  ans[[2]] <- fitmaturity(grp2,sams,length=length,mat10=mat10)
  ans[[3]] <- fitmaturity(c(grp1,grp2),sams,length=length,mat10=mat10)
  xmin <- getmin(c(ans[[1]]$plotdat[,"lens"],ans[[2]]$plotdat[,"lens"]))
  xmax <- getmax(c(ans[[1]]$plotdat[,"lens"],ans[[2]]$plotdat[,"lens"]))
  parset()
  plot(ans[[1]]$plotdat[,"lens"],ans[[1]]$plotdat[,"propm"],type="p",pch=16,
       cex=1.0,xlim=c(xmin,xmax),xlab="shell length",yaxs="r",
       ylab=paste0("Proportion Mature"),panel.first=grid())
  points(ans[[2]]$plotdat[,"lens"],ans[[2]]$plotdat[,"propm"],pch=1,
         col="blue",cex=1.75)
  lines(ans[[1]]$plotdat[,"lens"],ans[[1]]$plotdat[,"predm"],lwd=2,col=1)
  lines(ans[[2]]$plotdat[,"lens"],ans[[2]]$plotdat[,"predm"],lwd=2,col="blue")
  points(ans[[3]]$plotdat[,"lens"],ans[[3]]$plotdat[,"propm"],pch=1,col=2,
         cex=1.5)
  lines(ans[[3]]$plotdat[,"lens"],ans[[3]]$plotdat[,"predm"],lwd=2,col=2)
  lab1 <- makelabel("",grp1)
  lab2 <- makelabel("",grp2)
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
#' @param mat10 column name for the maturity variable in samdat mature=1,
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
#' print("wait on internal data sets")
#' # site1=grp;site2=grp2; sitecol="site";samdat=samb;length="length";mature="mature"
testgroup <- function(grp,samdat,sitecol="site",length="length",mat10="mature") { 
  picks <- which(samdat[,sitecol] %in% grp)
  sams <- samdat[picks,]
  samsf <- sams
  samsf[,sitecol] <- factor(samsf[,sitecol])
  model <- glm(formula = samsf[,mat10] ~ samsf[,length] + samsf[,sitecol], 
               family = binomial)
  summ <- summary(model)
  return(summ=summ)
} # end of testgroup
