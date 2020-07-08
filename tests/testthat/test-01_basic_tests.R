context("Does psupertime work?")
# devtools::document(); devtools::test()

################
# set up
################

suppressPackageStartupMessages({
  library('psupertime')
  library('SingleCellExperiment')
  })
seed <- as.numeric(format(Sys.time(), "%s"))
set.seed(seed)

# load the data
data(acinar_hvg_sce)

# run psupertime
y           = acinar_hvg_sce$donor_age


################
# tests
################

test_that("do standard calls to psupertime work?", {
    # do they work ok?
    expect_is(psupertime(acinar_hvg_sce, y, sel_genes='all'), 'list')
    expect_is(psupertime(acinar_hvg_sce, y, sel_genes='hvg'), 'list')
    expect_is(psupertime(acinar_hvg_sce, y, sel_genes=list(bio_cutoff=0, hvg_cutoff=1)), 'list')
})

test_that("do plotting functions work?", {
    # run psupertime
    psuper_obj  = psupertime(acinar_hvg_sce, y, sel_genes='all')

    # check each plotting function
    expect_is(plot_identified_gene_coefficients(psuper_obj), 'ggplot') 
    expect_is(plot_identified_genes_over_psupertime(psuper_obj), 'ggplot') 
    expect_is(plot_labels_over_psupertime(psuper_obj), 'ggplot') 
    expect_is(plot_predictions_against_classes(psuper_obj), 'ggplot') 
    expect_is(plot_train_results(psuper_obj), 'ggplot') 
})
