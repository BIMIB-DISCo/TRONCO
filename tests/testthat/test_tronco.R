
data(test_dataset)
data(test_dataset_no_hypos)
capri = tronco.capri(test_dataset, nboot = 1, silent = TRUE)
caprese = tronco.caprese(test_dataset_no_hypos, silent = TRUE)
chowliu = tronco.mst.chowliu(test_dataset_no_hypos, nboot = 1, silent = TRUE)
edmonds = tronco.mst.edmonds(test_dataset_no_hypos, nboot = 1, silent = TRUE)
prim = tronco.mst.prim(test_dataset_no_hypos, nboot = 1, silent = TRUE)

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
    expect_output(tronco.capri(test_dataset,
        nboot = 1,
        regularization = c('bic', 'aic', 'loglik')),
    "CAPRI")
    expect_output(tronco.capri(test_dataset,
        nboot = 1,
        regularization = c('bic'),
        command = 'tabu',
        do.boot = FALSE),
    "CAPRI")
    expect_output(tronco.capri(test_dataset,
        nboot = 1,
        boot.seed = 1),
    "CAPRI")
    expect_error(tronco.capri(NULL))
    expect_error(tronco.capri(test_dataset,
        command = 'banana'))
    expect_error(tronco.capri(test_dataset,
        pvalue = -1))
    expect_error(tronco.capri(test_dataset,
        pvalue = 2))
    expect_error(tronco.capri(test_dataset,
        regularization = 'banana'))
})

context("EDMONDS")

test_that("tronco edmonds is working", {
    expect_output(tronco.mst.edmonds(test_dataset_no_hypos,
        nboot = 1,
        regularization = c('no_reg', 'loglik', 'bic', 'aic')),
    "EDMONDS")
    expect_output(tronco.mst.edmonds(test_dataset_no_hypos,
        nboot = 1,
        boot.seed = 1),
    "EDMONDS")
    expect_warning(tronco.mst.edmonds(test_dataset,
        nboot = 1))
    expect_error(tronco.mst.edmonds(NULL))
    expect_error(tronco.mst.edmonds(test_dataset_no_hypos,
        command = 'banana'))
    expect_error(tronco.mst.edmonds(test_dataset_no_hypos,
        pvalue = -1))
    expect_error(tronco.mst.edmonds(test_dataset_no_hypos,
        pvalue = 2))
    expect_error(tronco.mst.edmonds(test_dataset_no_hypos,
        regularization = 'banana'))
})

context("CHOW LIU")

test_that("tronco chow liu is working", {
    expect_output(tronco.mst.chowliu(test_dataset_no_hypos,
        nboot = 1,
        regularization = c('loglik', 'bic', 'aic')),
    "Chow")
    expect_output(tronco.mst.chowliu(test_dataset_no_hypos,
        nboot = 1,
        boot.seed = 1),
    "Chow")
    expect_warning(tronco.mst.chowliu(test_dataset,
        nboot = 1))
    expect_error(tronco.mst.chowliu(NULL))
    expect_error(tronco.mst.chowliu(test_dataset_no_hypos,
        command = 'banana'))
    expect_error(tronco.mst.chowliu(test_dataset_no_hypos,
        pvalue = -1))
    expect_error(tronco.mst.chowliu(test_dataset_no_hypos,
        pvalue = 2))
    expect_error(tronco.mst.chowliu(test_dataset_no_hypos,
        regularization = 'banana'))
})

context("PRIM")

test_that("tronco prim is working", {
    expect_output(tronco.mst.prim(test_dataset_no_hypos,
        nboot = 1,
        regularization = c('no_reg', 'loglik', 'bic', 'aic')),
    "PRIM")
    expect_output(tronco.mst.prim(test_dataset_no_hypos,
        nboot = 1,
        boot.seed = 1),
    "PRIM")
    expect_warning(tronco.mst.prim(test_dataset,
        nboot = 1))
    expect_error(tronco.mst.prim(NULL))
    expect_error(tronco.mst.prim(test_dataset_no_hypos,
        command = 'banana'))
    expect_error(tronco.mst.prim(test_dataset_no_hypos,
        pvalue = -1))
    expect_error(tronco.mst.prim(test_dataset_no_hypos,
        pvalue = 2))
    expect_error(tronco.mst.prim(test_dataset_no_hypos,
        regularization = 'banana'))
})

context("BOOTSTRAP")

test_that("tronco bootstrap is working", {
    expect_output(tronco.bootstrap(capri, nboot = 1), 'non-parametric')
    expect_output(tronco.bootstrap(capri, nboot = 1, type = 'statistical'), 'statistical')
    expect_output(tronco.bootstrap(caprese, nboot = 1), 'non-parametric')
    expect_output(tronco.bootstrap(prim, nboot = 1), 'non-parametric')
    expect_output(tronco.bootstrap(prim, nboot = 1, type = 'statistical'), 'statistical')
    expect_output(tronco.bootstrap(edmonds, nboot = 1), 'non-parametric')
    expect_output(tronco.bootstrap(edmonds, nboot = 1, type = 'statistical'), 'statistical')
    expect_output(tronco.bootstrap(chowliu, nboot = 1), 'non-parametric')
    expect_output(tronco.bootstrap(chowliu, nboot = 1, type = 'statistical'), 'statistical')
})


