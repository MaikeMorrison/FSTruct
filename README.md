
<!-- README.md is generated from README.Rmd. Please edit that file -->

# FSTruct – Measure variability in population structure estimates

<!-- badges: start -->
<!-- badges: end -->

The FSTruct R package provides an Fst-based tool to compare the
variability of Q matrices (matrices with rows of ancestry coefficients,
the default output of population structure inference programs such a
STRUCTURE and ADMIXTURE).

The package includes four functions:

-   `Q_stat` calculates Fst/Fst^max (a normalized measure of
    variability) for a Q matrix.
-   `Q_bootstrap` generates bootstrap replicates of one or more Q
    matrices along with associated statistics, including Fst/Fst^max, as
    well as statistical tests to compare bootstrap distributions of
    Fst/Fst^max .
-   `Q_plot` plots Q matrices using `ggplot2`.
-   `Q_simulate` generates one or more Q matrices using the Dirichlet
    distribution.

This package accompanies the paper “FSTruct: an FST-based tool for
measuring ancestry variation in inference of population structure” by
Maike Morrison, Nicolas Alcala, and Noah Rosenberg. [You can access the
paper in *Molecular Ecology Resources* at this
link.](https://doi.org/10.1111/1755-0998.13647) Q matrix files for data
examples that appeared in this paper are available [on this GitHub
repository, as
well.](https://github.com/MaikeMorrison/FSTruct/tree/master/Q%20matrix%20files%20for%20examples%20in%20Morrison%20et%20al.%20(2022))

## Installation

<!-- You can install the released version of FSTruct from [CRAN](https://CRAN.R-project.org) with: -->
<!-- ``` r -->
<!-- install.packages("FSTruct") -->
<!-- ``` -->

You can install this package from GitHub with the following R code:

``` r
# install.packages("devtools")
devtools::install_github("MaikeMorrison/FSTruct")
```

**If you have not used R before:** [Install
R](https://www.r-project.org/) (not optional) and
[RStudio](https://www.rstudio.com/products/rstudio/) (optional). In the
R or RStudio console, paste `install.packages("devtools")` and hit enter
to run. Next run
`devtools::install_github("MaikeMorrison/FSTruct", dependenceies = TRUE)`
to install FSTruct. Then see the [introductory
vignette](https://maikemorrison.com/files/fstruct-intro) for details on
how to use this package.

## Usage

For a step-by-step guide to this package, read the [introductory
vignette](https://maikemorrison.com/files/fstruct-intro).

For details on the functions and their usage, read the [FSTruct
manual](https://maikemorrison.com/files/fstruct-manual.pdf).

<!-- ## Example -->
<!-- This is a basic example which shows you how to solve a common problem: -->
<!-- ```{r example} -->
<!-- library(fstruct) -->
<!-- ## basic example code -->
<!-- ``` -->
<!-- What is special about using `README.Rmd` instead of just `README.md`? You can include R chunks like so: -->
<!-- ```{r cars} -->
<!-- summary(cars) -->
<!-- ``` -->
<!-- You'll still need to render `README.Rmd` regularly, to keep `README.md` up-to-date. `devtools::build_readme()` is handy for this. You could also use GitHub Actions to re-render `README.Rmd` every time you push. An example workflow can be found here: <https://github.com/r-lib/actions/tree/master/examples>. -->
<!-- You can also embed plots, for example: -->
<!-- ```{r pressure, echo = FALSE} -->
<!-- plot(pressure) -->
<!-- ``` -->
<!-- In that case, don't forget to commit and push the resulting figure files, so they display on GitHub and CRAN. -->
