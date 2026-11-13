# magenta

[![Project Status: WIP – Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.](https://www.repostatus.org/badges/latest/wip.svg)](https://www.repostatus.org/#wip)
[![R-CMD-check](https://github.com/OJWatson/magenta/actions/workflows/check-standard.yaml/badge.svg)](https://github.com/OJWatson/magenta/actions/workflows/check-standard.yaml)
[![pkgdown](https://github.com/OJWatson/magenta/actions/workflows/pkgdown.yaml/badge.svg)](https://github.com/OJWatson/magenta/actions/workflows/pkgdown.yaml)


## Motivation

`magenta` is an individual-based simulation model of malaria epidemiology and parasite genetics, which was designed to extends the Imperial malaria model by tracking the infection history of individuals. With this addition, genetic characteristics of the parasite can be assessed for looking at both neutral genetic variation as well as loci under selection. 

The model is written in C++/Rcpp and interfaced with R, incorporating both the deterministic version of the Imperial College Transmission model developed by Joel Hellewell and Hannah Slater, which helps initialise the model near equilbrium. In addition it is used to quickly the mosquito feeding behaviour, which is passed to `magenta` in order to increase its speed.    

## Installation

`magenta` is installed from github with:


``` r
#install.packages("devtools")
devtools::install_github("OJWatson/magenta")
```

Note that on Windows, a toolkit (Rtools) needs to be installed separately for devtools to work. 

Once installed, the package can be loaded using:


``` r
# Load the package
library(magenta)
```

Documentation
-------------

*magenta* is documented on a [dedicated
 website](https://ojwatson.github.io/magenta).

This includes the following vignettes:

- [**Introduction to magenta**](https://ojwatson.github.io/magenta/articles/Introduction.html):
  general introduction to the design of `magenta` and how to use it
- [**Model output, loggers and saved state**](https://ojwatson.github.io/magenta/articles/Loggers.html):
  what a simulation returns, and how to save and restart runs
- [**Schedules, parameters and reproducibility**](https://ojwatson.github.io/magenta/articles/Schedule.html):
  the per-year schedule a run follows, the parameter lists and their
  validators, the result object, and what a seed reproduces
- [**Simulating antimalarial resistance**](https://ojwatson.github.io/magenta/articles/Resistance.html):
  resistance loci, drug definitions, prophylaxis, mutation, importation and fitness costs
- [**Drug deployment strategies**](https://ojwatson.github.io/magenta/articles/DrugStrategies.html):
  multiple first line therapies, drug cycling and non-malarial fevers
- [**Adding a new drug**](https://ojwatson.github.io/magenta/articles/NewDrugs.html):
  defining therapies the package does not ship with, worked through with
  single low-dose primaquine
- [**Checking magenta against the deterministic model**](https://ojwatson.github.io/magenta/articles/Deterministic.html):
  the mosquito rates under ITN/IRS and the EIR to prevalence and incidence
  relationships compared with the deterministic model the package carries

Vignettes are pre-evaluated. The editable sources are in `vignettes_src/` and are
built into `vignettes/` with `Rscript vignettes_src/build_vignettes.R`, so that building
the package never re-runs a simulation.

Contributors
------------
- [OJ Watson](https://github.com/ojwatson)
- [Bob Verity](https://github.com/bobverity)
- [Joel Hellewell](https://github.com/jhellewell14)
- [Hannah Slater](https://github.com/hannahslater)

**Maintainer:** OJ Watson (o.watson15@imperial.ac.uk)
