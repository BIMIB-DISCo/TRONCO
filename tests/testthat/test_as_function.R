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

context("as.genotypes")

test_that("as.genotypes returns a genotypes matrix", {
    expect_equal(length(as.genotypes(muts)), 39)
    expect_equal(as.genotypes(NULL), NULL)
})

context("as.samples")

test_that("as.samples returns a list of samples", {
    expect_equal(as.samples(muts), 
                 unique(as.character(maf$Tumor_Sample_Barcode)))
    expect_equal(as.samples(NULL), NULL)
})

context("as.genes")

test_that("as.genes returns a list of samples", {
    expect_equal(as.genes(muts), unique(as.character(maf$Hugo_Symbol)))
    expect_equal(as.genes(hypo), unique(as.character(maf$Hugo_Symbol)))
    expect_equal(as.genes(NULL), NULL)
    expect_error(as.genes(hypo, types = 'Pattern'))
    expect_equal(length(as.genes(hypo)), 13)
})

context("as.events")

test_that("as.events returns a list of samples", {
    expect_equal(length(as.events(muts)), 26)
    expect_equal(length(as.events(hypo, types = 'Mutation')), 26)
    expect_equal(as.events(hypo, keysToNames = TRUE)[1],
                 "Mutation A2BP1")
    expect_equal(as.events(NULL), NULL)
})

context("as.stages")

test_that("as.stages returns a list of stages", {
    expect_equal(length(as.stages(muts_stages)), 1)
    expect_equal(as.stages(muts), NA)
    expect_equal(as.stages(NULL), NA)
})

context("as.types")

test_that("as.types returns a list of types", {
    expect_equal(as.types(muts), 'Mutation')
    expect_equal(as.types(NULL), NULL)
})

context("as.colors")

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

context("as.gene")

test_that("as.gene returns a list of types", {
    expect_equal(as.gene(muts, genes='A2BP1')[1,1], 1)
})

context("as.alteration")

test_that("as.alteration returns a compliant dataset", {
    expect_equal(unique(as.alterations(gistic, silent = TRUE)$annotations[,'type']),
                 "Alteration")
    expect_error(as.alterations(hypo))
    expect_error(as.alterations(gistic_model))
})

context("as.patterns")

test_that("as.patterns return a list of patterns", {
    expect_null(as.patterns(muts))
    expect_equal(as.patterns(hypo),
                 "test")
})

context("as.hypotheses")

test_that("as.hypotheses return a list of hypotheses", {
    expect_null(as.hypotheses(gistic))
    expect_equal(length(as.hypotheses(hypo)), 88)
    expect_equal(length(as.hypotheses(hypo, cause = 'A2BP1')), 4)
    expect_equal(length(as.hypotheses(hypo, effect = 'A2BP1')), 4)
    expect_error(as.hypotheses(hypo, cause = 'X'))
    expect_error(as.hypotheses(hypo, effect = 'X'))
})

context("as.events.in.patterns")

test_that("as.events.in.patterns return the right amount of genes", {
    expect_equal(length(as.events.in.patterns(hypo)), 4)
    expect_equal(length(as.events.in.patterns(hypo, pattern = 'test')), 4)
    expect_error(as.events.in.patterns(hypo, pattern = 'X'))
})

context("as.genes.in.patterns")

test_that("as.genes.in.patterns return the right amount of genes", {
    expect_equal(length(as.genes.in.patterns(hypo)), 2)
})

context("as.types.in.patterns")

test_that("as.types.in.patterns return the right amount of genes", {
    expect_equal(length(as.types.in.patterns(hypo)), 1)
})

context("as.events.in.sample")

test_that("as.events.in.sample return the right amount of genes", {
    expect_equal(length(as.events.in.sample(muts, as.samples(muts)[1])), 8)
})

context("as.models")

test_that("as.models return a model", {
    expect_equal(length(as.models(gistic_model)), 1)
    expect_error(as.models(gistic_model, models = 'X'))
})

context("as.description")

test_that("as.description return a description", {
    expect_equal(as.description(hypo), "")
    desc = annotate.description(hypo, 'desc')
    expect_equal(as.description(desc), "desc")
})

context("as.pathway")

test_that("as.pathway return compliant dataset", {
    expect_output(as.pathway(muts,
                      pathway.name = 'test',
                      pathway.genes = 'APC'),
                "succesfully.")
    expect_output(as.pathway(muts,
                      pathway.name = 'test',
                      pathway.genes = 'APC',
                      aggregate.pathway = FALSE),
                "succesfully.")
    expect_output(as.pathway(muts_stages,
                      pathway.name = 'test',
                      pathway.genes = 'APC'),
                "succesfully.")
})

context("as.adj.matrix")

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

context("as.duplicates")

test_that("as.duplicates return false", {
    expect_false(has.duplicates(hypo))
})

context("as.parents.pos")

test_that("as.parents.pos return the right amount of models", {
    expect_equal(length(as.parents.pos(gistic_model_capri)), 2)
    expect_error(as.parents.pos(gistic_model_capri, models = 'banana'))
})

context("duplicates")

test_that("duplicates return 0", {
    expect_equal(length(duplicates(hypo)), 0)
})

context("npatterns")

test_that("npatterns return the right amount of pattern", {
    expect_equal(length(duplicates(hypo)), 0)
})

context("view")

test_that("view return output", {
    expect_output(view(gistic_model_capri), regexp=NULL)
})

context("order.frequency")

test_that("order.frequency reorder the genotypes", {
    expect_equal(length(order.frequency(gistic_model_capri)), 9)
})

context("enforce.string")

test_that("enforce.string return a string genotypes", {
    expect_equal(length(enforce.string(gistic_model_capri)), 9)
})
