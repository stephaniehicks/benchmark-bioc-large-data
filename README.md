
# Benchmarking large-scale data using `scry`

This repository is to benchmark the `nullResiduals()` function in `scry` [currently in the `hdf5` branch](https://github.com/kstreet13/scry/tree/hdf5). 
The `nullResiduals()` function performs a transformation that when combined with  PCA approximates GLM-PCA [implemented in the `glmpca`](https://github.com/willtownes/glmpca) R package. 

We compare the performance (time and memory) using:  

- data stored in-memory (RAM)
- data stored on-disk (using the `DelayedArray` framework and `HDF5Array` objects)

We consider two model types for calculating residuals: 

- `binomial` = (default) closest approximation to multinomial
- `poisson` =  may be faster to compute and often is very similar to `binomial`

We consider two types of residuals: 

- deviance residuals
- Pearson residuals

## Data 

We consider the following datasets: 

- [`TENxBrainData`](https://bioconductor.org/packages/release/data/experiment/html/TENxBrainData.html)


## Installation 

The [glmpca package](https://CRAN.R-project.org/package=glmpca) is available from CRAN and the [`scry` package](https://bioconductor.org/packages/scry) is available from Bioconductor. To install the stable releases (recommended):

```r
install.packages("glmpca")
BiocManager::install("scry")
```

To install the development versions (needed for this benchmark):

```r
remotes::install_github("willtownes/glmpca")
remotes::install_github("kstreet13/scry@hdf5") # install HDF5 branch
```

## Authors

- [Stephanie Hicks](https://www.stephaniehicks.com), https://github.com/stephaniehicks, Johns Hopkins Bloomberg School of Public Health, USA


## Issues and bug reports

Please use https://github.com/stephaniehicks/benchmark-bioc-large-data/issues to submit issues, bug reports, and comments.
