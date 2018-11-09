# Random correlation matrices with Metropolis-Hastings

## Contents

- `time_experiment.R`: script that executes related state-of-the-art algorithms
  to measure execution time.
- `plot_utils.R`: utility functions for plotting.
- `plot.R`: script that generates the plots describing the properties of the
  algorithm, as well as the results of the time experiment.

## Requirements
These scripts rely on the development branch `rchol` of the
`R` package `gmat`, which can be installed with the `devtools`
package from CRAN as follows

```R
# install.packages("devtools")
devtools::install_github("irenecrsn/gmat", ref = "rchol")
```
More information on the `gmat` package can be found on [its
repository](https://github.com/irenecrsn/gmat).

## Usage

```bash
	Rscript time_experiment.R
	Rscript plot.R
```
The time experiment can be computationally intensive.

