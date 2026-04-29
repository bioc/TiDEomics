
<!-- README.md is generated from README.Rmd. Please edit that file -->

# TiDEomics: Time-course Differential Expression analysis of omics data

<!-- badges: start -->

[![GitHub
issues](https://img.shields.io/github/issues/hte123/TiDEomics)](https://github.com/hte123/TiDEomics/issues)
[![GitHub
pulls](https://img.shields.io/github/issues-pr/hte123/TiDEomics)](https://github.com/hte123/TiDEomics/pulls)
[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- [![Bioc release status](http://www.bioconductor.org/shields/build/release/bioc/TiDEomics.svg)](https://bioconductor.org/checkResults/release/bioc-LATEST/TiDEomics) -->
<!-- [![Bioc devel status](http://www.bioconductor.org/shields/build/devel/bioc/TiDEomics.svg)](https://bioconductor.org/checkResults/devel/bioc-LATEST/TiDEomics) -->
<!-- [![Bioc downloads rank](https://bioconductor.org/shields/downloads/release/TiDEomics.svg)](http://bioconductor.org/packages/stats/bioc/TiDEomics/) -->
<!-- [![Bioc support](https://bioconductor.org/shields/posts/TiDEomics.svg)](https://support.bioconductor.org/tag/TiDEomics) -->
<!-- [![Bioc history](https://bioconductor.org/shields/years-in-bioc/TiDEomics.svg)](https://bioconductor.org/packages/release/bioc/html/TiDEomics.html#since) -->
<!-- [![Bioc last commit](https://bioconductor.org/shields/lastcommit/devel/bioc/TiDEomics.svg)](http://bioconductor.org/checkResults/devel/bioc-LATEST/TiDEomics/) -->
<!-- [![Bioc dependencies](https://bioconductor.org/shields/dependencies/release/TiDEomics.svg)](https://bioconductor.org/packages/release/bioc/html/TiDEomics.html#since) -->
<!-- [![check-bioc](https://github.com/hte123/TiDEomics/actions/workflows/check-bioc.yml/badge.svg)](https://github.com/hte123/TiDEomics/actions/workflows/check-bioc.yml) -->
<!-- [![Codecov test coverage](https://codecov.io/gh/hte123/TiDEomics/graph/badge.svg)](https://app.codecov.io/gh/hte123/TiDEomics) -->
<!-- badges: end -->

`TiDEomics` is designed to streamline **Ti**me-course **D**ifferential
**E**xpression analysis of **omics** data with **multiple experimental
groups / conditions** (for example, different cell lines or different
treatments sampled at several time points). The package’s main goals
are:

- Compare **multiple time courses** side-by-side (e.g. multiple mutants
  vs control).
- Identify differential expression by **time**, **group**, and **both**
  (time × group interactions).
- Provide utilities for QC, breakpoint detection (segmented regression),
  variance decomposition, and module detection (WGCNA) that are tailored
  for time-course multi-sample data.

## Installation

<!-- Get the latest stable `R` release from [CRAN](http://cran.r-project.org/). Then install `TiDEomics` from [Bioconductor](http://bioconductor.org/) using the following code: -->

<!-- ```{r 'install', eval = FALSE} -->

<!-- if (!requireNamespace("BiocManager", quietly = TRUE)) { -->

<!--     install.packages("BiocManager") -->

<!-- } -->

<!-- # BiocManager::install("TiDEomics") -->

<!-- ``` -->

Install the development version from
[GitHub](https://github.com/hte123/TiDEomics) with:

``` r
if (!require("remotes", quietly = TRUE)) install.packages("remotes")

remotes::install_github("hte123/TiDEomics")
```

## Example

For detailed examples and explanations, please refer to the
[tutorial](articles/TiDEomics.html), applications and other package
documentation.

## Citation

Below is the citation output from using `citation('TiDEomics')` in R.
Please run this yourself to check for any updates on how to cite
**TiDEomics**.

``` r
print(citation("TiDEomics"), bibtex = TRUE)
#> To cite package 'TiDEomics' in publications use:
#> 
#>   He T (2026). _TiDEomics: Time-course Differential Expression analysis
#>   of omics data_. R package version 0.99.0,
#>   <https://github.com/hte123/TiDEomics>.
#> 
#> A BibTeX entry for LaTeX users is
#> 
#>   @Manual{,
#>     title = {TiDEomics: Time-course Differential Expression analysis of omics data},
#>     author = {Tianen He},
#>     year = {2026},
#>     note = {R package version 0.99.0},
#>     url = {https://github.com/hte123/TiDEomics},
#>   }
```

Please note that the `TiDEomics` was only made possible thanks to many
other R and bioinformatics software authors, which are cited either in
the vignettes and/or the paper(s) describing this package.

## Code of Conduct

Please note that the `TiDEomics` project is released with a [Contributor
Code of Conduct](http://bioconductor.org/about/code-of-conduct/). By
contributing to this project, you agree to abide by its terms.

## Development tools

- Continuous code testing is possible thanks to [GitHub
  actions](https://www.tidyverse.org/blog/2020/04/usethis-1-6-0/)
  through *[usethis](https://CRAN.R-project.org/package=usethis)*,
  *[remotes](https://CRAN.R-project.org/package=remotes)*, and
  *[rcmdcheck](https://CRAN.R-project.org/package=rcmdcheck)* customized
  to use [Bioconductor’s docker
  containers](https://www.bioconductor.org/help/docker/) and
  *[BiocCheck](https://bioconductor.org/packages/3.23/BiocCheck)*.
- Code coverage assessment is possible thanks to
  [codecov](https://codecov.io/gh) and
  *[covr](https://CRAN.R-project.org/package=covr)*.
- The [documentation website](http://hte123.github.io/TiDEomics) is
  automatically updated thanks to
  *[pkgdown](https://CRAN.R-project.org/package=pkgdown)*.
- The code is styled automatically thanks to
  *[styler](https://CRAN.R-project.org/package=styler)*.
- The documentation is formatted thanks to
  *[devtools](https://CRAN.R-project.org/package=devtools)* and
  *[roxygen2](https://CRAN.R-project.org/package=roxygen2)*.

<!-- For more details, check the `dev` directory. -->

This package was developed using
*[biocthis](https://bioconductor.org/packages/3.23/biocthis)*.
