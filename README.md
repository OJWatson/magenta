
# magenta

[![Project Status: WIP – Initial development is in progress, but there
has not yet been a stable, usable release suitable for the
public.](https://www.repostatus.org/badges/latest/wip.svg)](https://www.repostatus.org/#wip)
[![Travis-CI Build
Status](https://travis-ci.org/ojwatson/magenta.png?branch=master)](https://travis-ci.org/ojwatson/magenta)

## Motivation

`magenta` is an individual-based simulation model of malaria
epidemiology and parasite genetics, which was designed to extends the
Imperial malaria model by tracking the infection history of individuals.
With this addition, genetic characteristics of the parasite can be
assessed for looking at both neutral genetic variation as well as loci
under selection.

The model is written in C++/Rcpp and interfaced with R, incorporating
both the deterministic version of the Imperial College Transmission
model developed by Joel Hellewell and Hannah Slater, which helps
initialise the model near equilibrium. In addition it is used to quickly
the mosquito feeding behaviour, which is passed to `magenta` in order to
increase its speed.

## Installation

You can also install `magenta` from github with:

``` r
#install.packages("devtools")
devtools::install_github("OJWatson/magenta")
```

Note that on Windows, a toolkit (Rtools) needs to be installed
separately for devtools to work.

Once installed, the package can be loaded using:

``` r
# Load the package
library(magenta)
```

## Documentation

*magenta* is documented on a [dedicated
website](https://ojwatson.github.io/magenta).

This includes the following vignettes:

  - **`Introduction`**: general introduction to the design of `magenta`
    and how to use it
  - **`Loggers`**: brief overview of the logging information returned by
    `magenta` (TO DO)
  - **`Resistance`**: (under development) description of how to simulate
    resistance

## Contributors

  - [OJ Watson](https://github.com/ojwatson)
  - [Bob Verity](https://github.com/bobverity)
  - [Joel Hellewell](https://github.com/jhellewell14)
  - [Hannah Slater](https://github.com/hannahslater)

**Maintainer:** OJ Watson (<oj.watson@hotmail.co.uk>)
