# Feature Extraction for HDLSS Labelled Data

## Contents

- [Overview](#overview)
- [Repo Contents](#repo-contents)
- [System Requirements](#system-requirements)
- [Installation Guide](#installation-guide)
- [Demo](#demo)
- [Results](#results)
- [License](./LICENSE)
- [Issues](https://github.com/ebridge2/fselect/issues)

# Overview

In modern scientific discovery, it is becoming increasingly critical to investigate data in the HDLSS (high dimensionality; low sample size) setting. In this package, we provide several methosd for feature extraction, reference simulations, and perform a comprehensive comparison.

# Repo Contents

- [R](./R): `R` package code.
- [docs](./docs): package documentation.
- [man](./man): package manual for help in R session.
- [tests](./tests): `R` unit tests written using the `testthat` package.
- [vignettes](./vignettes): `R` vignettes for R session html help pages.


# System Requirements

## Hardware Requirements

The `fselect` package requires only a standard computer with enough RAM to support the operations defined by a user. For minimal performance, this will be a computer with about 2 GB of RAM. For optimal performance, we recommend a computer with the following specs:

RAM: 16+ GB  
CPU: 4+ cores, 3.3+ GHz/core

The runtimes below are generated using a computer with the recommended specs (16 GB RAM, 4 cores@3.3 GHz) and internet of speed 25 Mbps.

## Software Requirements

### OS Requirements

This package is supported for *Linux* operating systems. The package has been tested on the following systems:

Linux: Ubuntu 16.04  
Mac OSX:  
Windows:  

Before setting up the `fselect` package, users should have `R` version 3.4.0 or higher, and several packages set up from CRAN.

#### Installing R version 3.4.2 on Ubuntu 16.04

the latest version of R can be installed by adding the latest repository to `apt`:

```
sudo echo "deb http://cran.rstudio.com/bin/linux/ubuntu xenial/" | sudo tee -a /etc/apt/sources.list
gpg --keyserver keyserver.ubuntu.com --recv-key E084DAB9
gpg -a --export E084DAB9 | sudo apt-key add -
sudo apt-get update
sudo apt-get install r-base r-base-dev
```

which should install in about 20 seconds.

#### Package dependencies

Users should install the following packages prior to installing `fselect`, from an `R` terminal:

```
install.packages(c('ggplot2', 'abind', 'irlba', 'knitr', 'rmarkdown', 'latex2exp', 'MASS'))
```

which will install in about 30 seconds on a recommended machine.

#### Package Versions

The `fselect` package functions with all packages in their latest versions as they appear on `CRAN` on December 13, 2017. Users can check [CRAN snapshot](https://mran.microsoft.com/timemachine/) for details. The versions of software are, specifically:
```
abind_1.4-5
latex2exp_0.4.0
ggplot2_2.2.1
irlba_2.3.1
Matrix_1.2-3
MASS_7.3-47
```

If you are having an issue that you believe to be tied to software versioning issues, please drop us an [Issue](https://github.com/neurodata/mgc/issues). 

# Installation Guide

From an `R` session, type:

```
require(devtools)
install_github('neurodata/fselect', build_vignettes=TRUE, force=TRUE)  # install fselect with the vignettes
require(fselect)
vignette("lol", package="fselect")  # view one of the basic vignettes
```

The package should take approximately 15 seconds to install with vignettes on a recommended computer. 

# Demo

For interactive demos of the functions, please check out the vignettes built into the package. They can be accessed as follows:

```
require(fselect)
vignette('lol')
vignette('pca')
vignette('cpca')
vignette('lrcca')
vignette('lda')
vignette('xval')
vignette('simulations')
```

# Results

[MNIST](https://htmlpreview.github.io/?https://github.com/neurodata/fselect/blob/master/docs/mnist.html)
