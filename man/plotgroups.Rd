% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/maturitycalc.r
\name{plotgroups}
\alias{plotgroups}
\title{plotgroups fits a maturity ~ length logistic regression to groups}
\usage{
plotgroups(
  groups,
  samdat,
  xmin = 0,
  xmax = 0,
  sitecol = "site",
  length = "length",
  mature = "mature",
  margin = c(0.4, 0.4, 0.05, 0.05)
)
}
\arguments{
\item{groups}{a list containing vectors of site numbers}

\item{samdat}{a data.frame of maturity data containing the required site(s)}

\item{xmin}{a generic x-axis minimum?, if 0 use local, default=0}

\item{xmax}{a generic x-axis maximum?, if 0 use local, default=0}

\item{sitecol}{column name of the variable for the site numbers. The default
column name = 'site'}

\item{length}{column name in samdat for the length variable}

\item{mature}{column name for the maturity variable in samdat mature=1,
immature = 0}

\item{margin}{what margin should be used in the plot? 
default=c(0.4,0.4,0.05,0.05)}
}
\value{
an invisible list containing the parameters and the fitmaturity 
    results. It also generates a plot.
}
\description{
plotgroups after using findgroups one can use plotgroups to 
    plot maturity curves to each group found. It uses fitgroups and outputs
    those results invisibly.
}
\examples{
args(plotgroups)
print("wait on internal data sets")
#  groups=group; samdat=samb;length="length";mature="mature";sitecol="site"
}
