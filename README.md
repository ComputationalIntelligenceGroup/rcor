# Random correlation matrices with Metropolis-Hastings

This repository contains the files for replicating the experiments described in
the paper

> Córdoba I., Varando G., Bielza C., Larrañaga P. A fast
Metropolis-Hastings method for generating random correlation matrices. Lecture Notes in
Computer Science (IDEAL 2018), vol 11314, pp. 117-124, 2018. 

## Contents

- `time_experiment.R`: script that executes related state-of-the-art algorithms
  to measure execution time.
- `plot_utils.R`: utility functions for plotting.
- `plot.R`: script that generates the plots describing the properties of the
  algorithm, as well as the results of the time experiment.

## R packages required
- These scripts rely on the development version of the
`R` package `gmat`, which can be installed with the `devtools`
package from CRAN as follows

```R
# install.packages("devtools")
devtools::install_github("irenecrsn/gmat")
```
Specifically, the function `gmat::chol_mh` implements the Metropolis-Hastings
algorithm described in the Córdoba et al. (2018), while `gmat::chol_polar`
implements the polar parametrization method by Pouhramadi and Wang (2015):

> Pourahmadi, M., Wang, X. Distribution of random correlation matrices:
Hyperspherical parameterization of the Cholesky factor, Statistics &
Probability Letters, 106:5-12, 2015.

This and more information on the `gmat` package can be found on [its
repository](https://github.com/irenecrsn/gmat).

- The `CRAN` package `clusterGeneration` is also necessary for executing the vine
and onion methods described by Lewandowski et al. (2009):

> Lewandowski, D., Kurowicka, D., Joe, H. Generating random correlation matrices based on vines and extended onion method,
Journal of Multivariate Analysis, 100:9, pp. 1989-2001, 2009.

- For plotting the graphics, it is necessary to have installed the `CRAN`
  packages `dplyr` and `ggplot2`.

## Usage

```bash
	Rscript time_experiment.R
	Rscript plot.R
```
The time experiment can be computationally intensive.

