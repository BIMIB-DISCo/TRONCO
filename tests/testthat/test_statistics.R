data(test_dataset)
data(test_dataset_no_hypos)

model = tronco.capri(test_dataset, nboot = 1, silent = TRUE)
eloss = tronco.kfold.eloss(model, k = 2, runs = 2, silent = TRUE)
prederr = tronco.kfold.prederr(model, k = 2, runs = 2, cores.ratio = 0, silent = TRUE)
posterr = tronco.kfold.posterr(model, k = 2, runs = 2, cores.ratio = 0, silent = TRUE)

context("eloss")

test_that("eloss", {
    expect_output(tronco.kfold.eloss(model, k = 2, runs = 2), 'loss')
    expect_equal(length(as.kfold.eloss(eloss)), 3) 
    expect_equal(length(as.confidence(eloss, conf = 'eloss')), 1)
})

context("prederr")

test_that("eloss", {
    expect_output(tronco.kfold.prederr(model, k = 2, runs = 2), 'prediction')
    expect_equal(length(as.kfold.prederr(prederr)), 2)
    expect_equal(length(as.kfold.prederr(prederr, table = TRUE)), 2)
    expect_equal(length(as.confidence(prederr, conf = 'prederr')), 1)
})


context("posterr")

test_that("posterr", {
    expect_output(tronco.kfold.posterr(model, k = 2, runs = 2), 'posterior')
    expect_equal(length(as.kfold.posterr(posterr)), 2)
    expect_equal(length(as.kfold.posterr(posterr, table = TRUE)), 2)
    expect_equal(length(as.confidence(posterr, conf = 'posterr')), 1)
})
