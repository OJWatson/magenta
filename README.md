# magenta

[![Project Status: WIP – Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.](https://www.repostatus.org/badges/latest/wip.svg)](https://www.repostatus.org/#wip)


## Motivation

`magenta` is an individual-based simulation model of malaria epidemiology and parasite genetics, which was designed to extends the Imperial malaria model by tracking the infection history of individuals. With this addition, genetic characteristics of the parasite can be assessed for looking at both neutral genetic variation as well as loci under selection. 

The model is written in C++/Rcpp and interfaced with R, incorporating both the deterministic version of the Imperial College Transmission model developed by Joel Hellewell and Hannah Slater, which helps initialise the model near equilbrium. In addition it is used to quickly the mosquito feeding behaviour, which is passed to `magenta` in order to increase its speed.    

## Installation

You can install the development version of `magenta` from GitHub with:


```r
#install.packages("devtools")
devtools::install_github("OJWatson/magenta")
```

Note that on Windows, a toolkit (Rtools) needs to be installed separately for devtools to work. 

Once installed, the package can be loaded using:


```r
# Load the package
library(magenta)
```

## Documentation

*magenta* is documented on a [dedicated website](https://ojwatson.github.io/magenta).

Key articles currently available are:

- **`Introduction`** – overview of the simulator and first steps for running experiments.
- **`Loggers`** – details on the logging outputs, update intervals, and how to persist model state.
- **`Resistance`** – guidance on configuring drug resistance scenarios and interpreting outputs.
- **`Drug components`** – deeper discussion of the drug parameterisation and helper utilities.
- **`Identity-by-descent (IBD)`** – explanation of the IBD mode and how to interpret IBD summaries.
- **`Advanced features`** – notes on modelling non-malarial fevers, vector adaptation, and sexual parasite stages.

Each vignette lives in `vignettes/` and can be browsed locally via `browseVignettes("magenta")` after installation.

## Contributors
- [OJ Watson](https://github.com/ojwatson)
- [Bob Verity](https://github.com/bobverity)
- [Joel Hellewell](https://github.com/jhellewell14)
- [Hannah Slater](https://github.com/hannahslater)

**Maintainer:** OJ Watson (oj.watson@imperial.ac.uk)
