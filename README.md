
<!-- README.md is generated from README.Rmd. Please edit that file -->

# fstruct – Measure variability in population structure estimates

<!-- badges: start -->

<!-- badges: end -->

The fstruct package provides a tool to compare variability between \(Q\)
matrices with rows of individual ancestry coefficients (the default
output of population structure inference programs such as STRUCTURE and
ADMIXTURE). The method employs the population differentiation statistic
\[\(f_{ST}\)\]{<https://en.wikipedia.org/wiki/Fixation_index>}.

The package includes three functions:

  - `structure.plot` plots \(Q\) matrices using `ggplot2`
  - `fst_stat` calculates \(F_{ST}/F_{ST}^\text{max}\) (a normalized
    measure of variability) for a \(Q\) matrix
  - `bootstrap` generates bootstrap replicates of one or more \(Q\)
    matrices along with associated statistics, including
    \(F_{ST}/F_{ST}^\text{max}\), as well as statistical tests to
    compare bootstrap distributions of \(F_{ST}/F_{ST}^\text{max}\) .

This package accompanies (insert paper citation here).

## Installation

<!-- You can install the released version of fstruct from [CRAN](https://CRAN.R-project.org) with: -->

<!-- ``` r -->

<!-- install.packages("fstruct") -->

<!-- ``` -->

You can install the development version from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("MaikeMorrison/fstruct")
```

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
