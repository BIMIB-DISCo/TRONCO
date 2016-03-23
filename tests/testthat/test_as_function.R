context("AS functions test")

data(maf)
muts = import.MAF(maf)
hypo = hypothesis.add(muts, 'test', OR('ABAT', 'ABCC3'))
no_hypo = delete.hypothesis(hypo, 'test')
data(gistic)
gistic = import.GISTIC(gistic)
gistic_model = tronco.caprese(gistic)
gistic_model_capri = tronco.capri(gistic, nboot = 2)

stages = list()
stages$stage = c('A', 'B', 'C')
names(stages$stage) = as.samples(muts)
stages = as.data.frame(stages)
muts_stages = annotate.stages(muts, stages=stages)

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
    expect_error(as.genes(hypo, types = 'Pattern'))
    expect_equal(length(as.genes(hypo)), 13)
})

test_that("as.events returns a list of samples", {
    data(as.events.test)
    expect_equal(as.events(muts), as.events.test)
    expect_equal(length(as.events(hypo, types = 'Mutation')), 26)
    expect_equal(as.events(hypo, keysToNames = TRUE)[1],
                 "Mutation A2BP1")
    expect_equal(as.events(NULL), NULL)
})

test_that("as.stages returns a list of stages", {
    data(as.stages.test)
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

test_that("as.events.in.patterns return the right amount of genes", {
    expect_equal(length(as.events.in.patterns(hypo)), 4)
    expect_equal(length(as.events.in.patterns(hypo, pattern = 'test')), 4)
    expect_error(as.events.in.patterns(hypo, pattern = 'X'))
})

test_that("as.genes.in.patterns return the right amount of genes", {
    expect_equal(length(as.genes.in.patterns(hypo)), 2)
})

test_that("as.types.in.patterns return the right amount of genes", {
    expect_equal(length(as.types.in.patterns(hypo)), 1)
})

test_that("as.events.in.sample return the right amount of genes", {
    expect_equal(length(as.events.in.sample(muts, as.samples(muts)[1])), 8)
})

test_that("as.models return a model", {
    expect_equal(length(as.models(gistic_model)), 1)
    expect_error(as.models(gistic_model, models = 'X'))
})

test_that("as.description return a description", {
    expect_equal(as.description(hypo), "")
    desc = annotate.description(hypo, 'desc')
    expect_equal(as.description(desc), "desc")
})

test_that("as.pathway return compliant dataset", {
    path = as.pathway(muts,
                      pathway.name = 'test',
                      pathway.genes = 'APC')
    expect_silent(is.compliant(path))
    path = as.pathway(muts,
                      pathway.name = 'test',
                      pathway.genes = 'APC',
                      aggregate.pathway = TRUE)
    expect_silent(is.compliant(path))
    path = as.pathway(muts_stages,
                      pathway.name = 'test',
                      pathway.genes = 'APC')
    expect_silent(is.compliant(path))
})

test_that("as.adj.matrix return a matrix", {
    expect_equal(length(as.adj.matrix(gistic_model_capri, 
                                      type = 'fit')), 2)
    expect_equal(length(as.adj.matrix(gistic_model_capri,
                                      type = 'pf')), 2)
    expect_error(as.adj.matrix(gistic_model_capri,
                               type = 'X'))
    expect_error(as.adj.matrix(gistic_model_capri,
                               models = 'X'))
})

test_that("as.duplicates return false", {
    expect_false(has.duplicates(hypo))
})

test_that("duplicates return 0", {
    expect_equal(length(duplicates(hypo)), 0)
})

test_that("npatterns return the right amount of pattern", {
    expect_equal(length(duplicates(hypo)), 0)
})


