data(maf)
muts = import.MAF(maf, silent = TRUE)
hypo = hypothesis.add(muts, 'test', OR('ABAT', 'ABCC3'))
no_hypo = delete.hypothesis(hypo, 'test')
data(crc_gistic)
gistic_err = import.GISTIC(crc_gistic, silent = TRUE)
gistic = delete.gene(gistic_err, 'APC')
gistic_model = tronco.caprese(gistic, silent = TRUE)
gistic_model_capri = tronco.capri(gistic, nboot = 2, silent = TRUE)

stages = list()
stages$stage = c('A', 'B', 'C')
names(stages$stage) = as.samples(muts)
stages = as.data.frame(stages)
muts_stages = annotate.stages(muts, stages=stages)

context("oncoprint")

test_that("oncoprint accept correct data", {
    expect_output(oncoprint(muts), "Oncoprint")
})

context("tronco.plot")

test_that("tronco.plot accept correct data", {
    expect_output(tronco.plot(gistic_model_capri), "Plotting")
})