
data(test_dataset)
data(test_dataset_no_hypos)

context("CAPRESE")

test_that("tronco caprese is working", {
    expect_output(tronco.caprese(test_dataset_no_hypos), "CAPRESE")
    expect_warning(tronco.caprese(test_dataset))
    expect_error(tronco.caprese(NULL))
    expect_error(tronco.caprese(test_dataset_no_hypos, lambda = -1))
    expect_error(tronco.caprese(test_dataset_no_hypos, lambda = 2))
})

context("CAPRI")

test_that("tronco capri is working", {
	expect_output(tronco.capri(test_dataset, nboot = 1), "CAPRI")
	expect_output(tronco.capri(test_dataset, nboot = 1, boot.seed = 1), "CAPRI")
	expect_error(tronco.capri(NULL))
	expect_error(tronco.capri(test_dataset, command = 'banana'))
	expect_error(tronco.capri(test_dataset, pvalue = -1))
    expect_error(tronco.capri(test_dataset, pvalue = 2))
    expect_error(tronco.capri(test_dataset, regularization = 'banana'))
})
