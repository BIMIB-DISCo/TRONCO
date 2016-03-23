data(maf)
muts = import.MAF(maf)
hypo = hypothesis.add(muts, 'test', OR('ABAT', 'ABCC3'))
no_hypo = delete.hypothesis(hypo, 'test')
data(gistic)
gistic = import.GISTIC(gistic)
gistic_model = tronco.caprese(gistic)


context("AS functions test")

test_that("as.genotypes returns a genotypes matrix", {
    data(as.genotypes.test)
    expect_equal(as.genotypes(muts), as.genotypes.test)
    expect_equal(as.genotypes(NULL), NULL)
})

test_that("as.samples returns a list of samples", {
    expect_equal(as.samples(muts), 
                 unique(as.character(maf$Tumor_Sample_Barcode)))
    expect_equal(as.samples(NULL), NULL)
})

test_that("as.genes returns a list of samples", {
    expect_equal(as.genes(muts), unique(as.character(maf$Hugo_Symbol)))
    expect_equal(as.genes(hypo), unique(as.character(maf$Hugo_Symbol)))
    expect_equal(as.genes(NULL), NULL)
})

test_that("as.events returns a list of samples", {
    data(as.events.test)
    expect_equal(as.events(muts), as.events.test)
    expect_equal(as.events(hypo, keysToName = TRUE)[1],
                 "Mutation A2BP1")
    expect_equal(as.events(NULL), NULL)
})

test_that("as.stages returns a list of stages", {
    data(as.stages.test)
    stages = list()
    stages$stage = c('A', 'B', 'C')
    names(stages$stage) = as.samples(muts)
    stages = as.data.frame(stages)
    muts_stages = annotate.stages(muts, stages=stages)
    expect_equal(as.stages(muts_stages), as.stages.test)
    expect_equal(as.stages(muts), NA)
    expect_equal(as.stages(NULL), NA)
})

test_that("as.types returns a list of types", {
    expect_equal(as.types(muts), 'Mutation')
    expect_equal(as.types(NULL), NULL)
})

test_that("as.colors returns a list of types", {
    areColors <- function(x) {
        sapply(x, function(X) {
            tryCatch(is.matrix(col2rgb(X)), 
                error = function(e) FALSE)
        })
    }
    expect_true(all(areColors(as.colors(muts))))
    expect_equal(as.colors(NULL), NULL)
})

test_that("as.gene returns a list of types", {
    expect_equal(as.gene(muts, genes='A2BP1')[1,1],
                 1)
})


test_that("as.alteration returns a compliant dataset", {
    expect_equal(unique(as.alterations(gistic)$annotations[,'type']),
                 "Alteration")
    expect_error(as.alterations(hypo))
    expect_error(as.alterations(gistic_model))
})

test_that("as.patterns return a list of patterns", {
    expect_null(as.patterns(muts))
    expect_equal(as.patterns(hypo),
                 "test")
})

test_that("as.hypotheses return a list of hypotheses", {
    expect_null(as.hypotheses(gistic))
    expect_equal(length(as.hypotheses(hypo)), 88)
    expect_equal(length(as.hypotheses(hypo, cause = 'A2BP1')), 4)
    expect_equal(length(as.hypotheses(hypo, effect = 'A2BP1')), 4)
    expect_error(as.hypotheses(hypo, cause = 'X'))
    expect_error(as.hypotheses(hypo, effect = 'X'))
})

