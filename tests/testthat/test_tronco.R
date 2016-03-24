context("TRONCO functions test")

data(test_dataset)
data(test_dataset_no_hypos)

test_that("tronco caprese is working", {
    expect_output(tronco.caprese(test_dataset_no_hypos), "CAPRESE")
    expect_warning(tronco.caprese(test_dataset))
    expect_error(tronco.caprese(NULL))
    expect_error(tronco.caprese(test_dataset_no_hypos, lambda = -1))
    expect_error(tronco.caprese(test_dataset_no_hypos, lambda = 2))
})
