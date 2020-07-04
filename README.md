
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

  - 28/06/2020 Initial combination of the older packages invLogistic and
    maturity.

  - 04/07/2020 added extra functionality for plots and analyses.

Malcolm Haddon

Hobart, June 28, 2020

## 

Haddon, M., Mundy, C., and D. Tarbath (2008) Using an inverse-logistic
model to describe growth increments of blacklip abalone (*Haliotis
rubra*) in Tasmania. *Fishery Bulletin* **106**:58-71
