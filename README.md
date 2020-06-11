# Random correlation matrices with Metropolis-Hastings

This repository contains the files for replicating the experiments described in
the paper

> Córdoba I., Varando G., Bielza C., Larrañaga P. A fast
Metropolis-Hastings method for generating random correlation matrices. Lecture Notes in
Computer Science (IDEAL 2018), vol 11314, pp. 117-124, 2018. 

Note that the time results described in the above paper may vary since we now
use package `randcorr` for the polar parametrization, which is now computationally
competitive with the rest of the compared methods.

## Contents

- `time_experiment.R`: script that executes related state-of-the-art algorithms
  to measure execution time.
- `plot_utils.R`: utility functions for plotting.
- `plot.R`: script that generates the plots describing the properties of the
  algorithm, as well as the results of the time experiment.
- `plot_comp.R`: optional script for generating additional comparison results.

## R packages required

These scripts rely on the R packages `gmat`, `clusterGeneration` and
`randcorr` for executing the algorithms described in Córdoba et al. (2018). For
plotting the graphics, it is also necessary to have installed the packages
`dplyr` and `ggplot2`.

## Usage

```bash
Rscript time_experiment.R
Rscript plot.R
```

