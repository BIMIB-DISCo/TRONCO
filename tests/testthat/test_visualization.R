data(maf)
muts = import.MAF(maf, silent = TRUE)
stages = list()
stages$stage = c('A', 'B', 'C')
names(stages$stage) = as.samples(muts)
stages = as.data.frame(stages)
muts_stages = annotate.stages(muts, stages=stages)

data(test_model_kfold)

context("oncoprint")

test_that("oncoprint accept correct data", {
    expect_output(oncoprint(muts), "Oncoprint")
})

context("tronco.plot")

test_that("tronco.plot accept correct data", {
    expect_output(tronco.plot(test_model_kfold), "Plotting")
})