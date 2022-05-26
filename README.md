
<!-- README.md is generated from README.Rmd. Please edit that file -->

# biology

An R package to faciliate the fitting and exploration of growth models
to invertebrate tagging data. In particular there is emphasis on the
inverse logistic curve, but also functions relating to the dose-response
and Michaelis-Menton curves.

To install this development version from
[GitHub](https://www.github.com/) you can use:

``` r
if (!require(devtools)){install.packages("devtools")} 

devtools::install_github("https://github.com/haddonm/biology")
```

There are currently no vignettes, although each function has its own
help page.

If you were to type *methods(“plot”)* into the console you will find
both plot.bootIL and plot.IL in the list. There are help pages for
those, as well as for print.IL and summary.IL. To see the contents of
the objects output by, say, *fitIL*, you can use, for example,
*str(ans)*.

Some details of the package can be found using
*packageDescription(“biology”)*.

-   26/05/2022 0.0.0.500 Many changes, still under development.

-   02/05/2022 Still needs a vignette and an internal data set, but can
    now be used to conduct an analysis of heterogeniety of maturity in a
    maturity data.set. It can group and separate sites in a
    statistically valid manner. No EWN

-   28/04/2022 added further functions to assist with working with
    groups of sites.

-   26/04/2022 Added many functions assisting the fitting of maturity
    curves

-   28/06/2020 Initial combination of the older packages invLogistic and
    maturity.

-   04/07/2020 added extra functionality for plots and analyses.

-   12/07/2020 0.0.0.900 Modified plot.IL and plotmodelIL. Amended
    fitDR, but still have to develop code to use the dose-response
    curve. No Errors, Warnings, or Notes

-   07/08/2020 0.0.0.850 developed the use of alternative dose response
    curve, using two alternative functions for the standard deviation.
    The first was constant and the second was an inverse relationship
    with initial length.

-   09/08/2020 0.0.0.800 Implemented a general funciton for fitting
    growth curves to taggign growth data. It needs a function to define
    which growth fucntion should be fitted, and a function to describe
    how the variance around the predicted growth increments should
    change as the predicted growth increments change.

Malcolm Haddon

Hobart, June 28, 2020

## 

Haddon, M., Mundy, C., and D. Tarbath (2008) Using an inverse-logistic
model to describe growth increments of blacklip abalone (*Haliotis
rubra*) in Tasmania. *Fishery Bulletin* **106**:58-71
