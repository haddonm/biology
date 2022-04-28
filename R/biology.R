

#' @importFrom grDevices palette rgb
#' @importFrom graphics abline grid hist lines mtext par plot points
#' @importFrom graphics text title legend polygon
#' @importFrom stats dnorm nlm optim quantile sd coef binomial glm anova
#' @importFrom hutils getmax getmin makelabel
#' @importFrom hplot parset pickbound plotprep
NULL


#' @title biology functions for Invertebrate Maturity and Tagging Growth
#'
#' @description The biology package provides functions to facilitate
#'     the analysis of invertebrate maturity data and invertebrate
#'     tagging growth increment data. Emphasis is given to using the 'inverse 
#'     logistic', a tagging data based growth curve developed by Haddon et al. 
#'     (2008), although alternatives are also implemented. As the name suggests 
#'     it predicts approximately linear growth (constant growth increments) in 
#'     juvenile stages,followed by reducing increments as the initial size/length
#'     increases, until the growth increments become very small. One advantage 
#'     of this curve is that it avoids predicting negative growth increments. A 
#'     development version is available on GitHub at github.com/haddonm/biology.
#'     The analysis of maturity curves (maturity-at-size or maturity-at-age) is
#'     based around classical mature/immature binomial data and uses a GLM with
#'     binomial errors. Other factors, such as site, or year, can be included in
#'     the analysis to determine the extent of variation through location and 
#'     time.
#'
#' @references Haddon, M., Mundy, C., and D. Tarbath (2008) Using an
#'     inverse-logistic model to describe growth increments of blacklip
#'     abalone (\emph{Haliotis rubra}) in Tasmania. \emph{Fishery Bulletin}
#'     \bold{106}: 58-71.
#'
#' @section Plotting and printing functions:
#' \describe{
#'   \item{addnorm}{ adds a normal distribution to the output from hist}
#'   \item{addlnorm}{ adds a log-normal distribution to output from hist}
#'   \item{inthist}{ plots a histogram of integer values more precisely
#'      than hist.}
#'   \item{newplot}{ is a bare-bones setup routine to generate a plot
#'      in RStudio using a floating window. If you want to alter the
#'      default par settings then you can use either setplot() to get
#'      suitable syntax or, more simply, use parsyn() which only give
#'      a template for the par syntax}
#'   \item{parset}{ defines the par statement for a single plot}
#'   \item{parsyn}{ types the standard syntax for the par command to
#'      the console}
#'   \item{plot1}{ simplifies the plotting of two variables in a single
#'      plot}
#'   \item{plot2}{ sets up a plotting window for two plots}
#'   \item{plotprep}{ Sets up a window and the par values for a plot.
#'      it checks to see if a graphics device is open and opens a new
#'      one if not. This is simply a utility function to save typing
#'      the standard syntax. Some of the defaults can be changed.
#'      Typing the name without () will provide a template for
#'      modification. If 'windows' is called repeatedly this will
#'      generate a new active graphics device each time leaving the
#'      older ones inactive but present. For quick exploratory plots
#'      this behaviour is not wanted, hence the check if an active
#'      device exists already or not.}
#'   \item{printV}{ returns a vector cbinded to 1:length(invect),
#'      which effectively prints the numbers vertically}
#' }
#' @section Data sets:
#' \describe{
#'   \item{midg}{a tagging data-set from abalone taken from the Middle
#'       Ground in the Actaeons.}
#'   \item{tasab}{Abalone maturity data from two sites on the
#'       south-west of Tasmania.}
#' }
#' @docType package
#' @name biology
NULL



