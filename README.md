[![Build Status](https://travis-ci.org/wmacnair/psupertime.svg?branch=master)](https://travis-ci.org/wmacnair/psupertime)
<br />
<img src="inst/extdata/psuper_sticker.jpeg" width="200" align="right"/> 

# psupertime

:wave: Hello User! :wave:

`psupertime` is an R package which uses single cell RNAseq data, where the cells have labels following a known sequence (e.g. a time series), to identify a small number of genes which place cells in that known order. It can be used for discovery of relevant genes, for exploration of unlabelled data, and assessment of one dataset with respect to the labels known for another dataset.

Read the pre-print here:
https://www.biorxiv.org/content/10.1101/622001v1

## How to install / use

To use this development version of the package, run the following lines in R:
```R
remotes::install_github('wmacnair/psupertime', build = TRUE, build_opts = c("--no-resave-data", "--no-manual"))
library('psupertime')
```
(You may need to install the package `remotes`, with `install.packages('remotes')`. Installation took <90s on a Macbook Pro.)

This should load all of the code and relevant documentation. 

## Basic analyses

We have included a small dataset which allows you to use some of the basic functionality in `psupertime`. To do this, have a look at the vignettes:
```R
browseVignettes(package = 'psupertime')
```
`psupertime` is fast: running these analyses took a bit under 1 minute on a Macbook Pro. 
The vignette also describes some of the additional functionality you can use, and full details are given in the documentation, via ```?psupertime```.


## Replicating analyses in the manuscript

To keep this main package light, we have only included a small example dataset. To replicate the figures in the manuscript and provide additional datasets for user experimentation, we have also made a data package, `psupplementary`. If you would like to see in more detail what `psupertime` can do, please go [here](https://github.com/wmacnair/psupplementary).


## Development roadmap

At the moment, `psupertime` takes a `SingleCellExperiment`/`sce` object as an input, and returns a `psupertime` object. I want to make this a bit smoother for users, by integrating the outputs from `psupertime` into the row and column annotations of the `sce`. So for example, the gene coefficients would be stored in `rowData`, and the latent time estimates would be stored in `colData`. I'll also add some unit tests and submit to `Bioconductor`.

## Suggestions

Please add any issues or requests to the _Issues_ page. All feedback enthusiastically received.

Cheers

Will




### System requirements

`psupertime` requires R (>= 3.4.3), and the following dependencies: `ggplot2` 3.1.1, `data.table` 1.12.2, `glmnet` 2.0-16, `scales` 1.0.0, `stringr` 1.4.0, `scran` 1.10.2, `SingleCellExperiment` 1.4.1, `SummarizedExperiment` 1.12.0, `RColorBrewer` 1.1-2.

