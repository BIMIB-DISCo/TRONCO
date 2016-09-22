data(maf)
muts = import.MAF(maf, silent = TRUE)
stages = list()
stages$stage = c('A', 'B', 'C')
names(stages$stage) = as.samples(muts)
stages = as.data.frame(stages)
muts_stages = annotate.stages(muts, stages=stages)

data(test_model_kfold)
data(test_model)
data(test_dataset)

context("oncoprint")

test_that("oncoprint accept correct data", {
    expect_output(oncoprint(muts), "Oncoprint")
})

context("tronco.plot")

test_that("tronco.plot accept correct data", {
    expect_output(tronco.plot(test_model), "Plotting")
    expect_output(tronco.plot(test_model, samples.annotation = 'patient 1'), "sample")
    expect_output(tronco.plot(test_model, pf = TRUE), "Plotting")
    expect_output(tronco.plot(test_model, expand = FALSE), "Plotting")
    expect_output(tronco.plot(test_model, fontsize = 20, scale.nodes = 0.3), "Plotting")
    expect_output(tronco.plot(test_model_kfold, conf = c('tp', 'pr', 'hg', 'eloss', 'prederr', 'posterr', 'npb')), "confidence")
    expect_output(tronco.plot(test_model, export.igraph = TRUE), "Plotting")
    expect_error(tronco.plot(test_dataset))
    expect_error(tronco.plot(test_model, models = c('a', 'b', 'c')))
    expect_error(tronco.plot(test_model, models = 'banana'))
    expect_error(tronco.plot(test_model, models = c('capri_bic', 'banana')))
    expect_error(tronco.plot(test_model, samples.annotation = c('s1', 's2'), pathways = TRUE))
    expect_error(tronco.plot(test_model, samples.annotation = c('s1', 's2')))
})
