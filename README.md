# psupertime

:wave: Hello User! :wave:

`psupertime` is an R package which uses single cell RNAseq data, where the cells have labels following a known sequence (e.g. a time series), to identify a small number of genes which place cells in that known order. It can be used for discovery of relevant genes, for exploration of unlabelled data, and assessment of one dataset with respect to the labels known for another dataset.


## How to install / use

To use this development version of the package, run the following lines in R:
```R
devtools::install_github('wmacnair/psupertime')
library('psupertime')
```

This should load all of the code and relevant documentation. 

## Basic analyses

We have included a small dataset which allows you to use some of the basic functionality in `psupertime`. To do this, have a look at the vignettes:
```R
browseVignettes(package = 'psupertime')
```

## Replicating analyses in the manuscript

To keep this main package light, we have only included a small example dataset. To replicate the figures in the manuscript and provide additional datasets for user experimentation, we have also made a data package, `psupplementary`. If you would like to see in more detail what `psupertime` can do, please go [here](https://github.com/wmacnair/psupplementary).


## Suggestions

Please add any issues or requests to the _Issues_ page. All feedback enthusiastically received.

Cheers

Will
