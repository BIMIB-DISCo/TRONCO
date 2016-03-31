data(test_dataset_no_hypos)
data(test_model)
data(test_dataset)

context("change.color")

test_that("change.color is working", {
    expect_equal(length(change.color(test_dataset, 'ins_del', 'red')), 5)
    expect_error(change.color(test_dataset_no_hypos, "XXX", "red"))
})

context("rename.type")

test_that("rename.type is working", {
    expect_equal(length(rename.type(test_dataset, 'ins_del', 'banana')), 5)
    expect_equal(length(rename.type(test_dataset, 'ins_del', 'ins_del')), 5)
    expect_equal(length(rename.type(test_dataset, 'ins_del', 'missense_point')), 5)
    expect_error(rename.type(test_dataset, 'Pattern', 'missense_point'))
    expect_error(rename.type(test_dataset, 'ins_del', 'Pattern'))
    expect_error(rename.type(test_dataset, 'banana', 'Pattern'))
})

context("rename.gene")

test_that("rename.gene is working", {
    expect_equal(length(rename.gene(test_dataset, 'TET2', 'banana')), 5)
    expect_error(rename.gene(test_dataset, 'XXX', 'banana'))
})

context("rename.gene")

test_that("rename.gene is working", {
    expect_equal(length(rename.gene(test_dataset, 'TET2', 'banana')), 5)
    expect_error(rename.gene(test_dataset, 'XXX', 'banana'))
})

context("delete.type")

test_that("delete.type is working", {
    expect_equal(length(delete.type(test_dataset_no_hypos, 'ins_del')), 5)
    expect_equal(length(delete.type(test_dataset, 'Pattern')), 5)
    expect_error(delete.type(test_model, 'ins_del'))
    expect_error(delete.type(test_dataset, 'ins_del'))
    expect_error(delete.type(test_dataset_no_hypos, 'banana'))
})

context("delete.gene")

test_that("delete.gene is working", {
    expect_equal(length(delete.gene(test_dataset_no_hypos, 'TET2')), 5)
    expect_error(delete.gene(test_model, 'TET2'))
    expect_error(delete.gene(test_dataset, 'EZH2'))
    expect_error(delete.gene(test_dataset_no_hypos, 'XXX'))
})

context("delete.event")

test_that("delete.event is working", {
    expect_equal(length(delete.event(test_dataset_no_hypos, 'TET2', 'ins_del')), 5)
    expect_error(delete.event(test_dataset, 'TET2', 'Pattern'))
    expect_error(delete.event(test_model, 'TET2', 'ins_del'))
    expect_error(delete.event(test_dataset, 'EZH2', 'ins_del'))
    expect_error(delete.event(test_dataset, 'XXX', 'Pattern'))
})

context("delete.hypothesis")

test_that("delete.hypothesis is working", {
    expect_equal(length(delete.hypothesis(test_dataset, event = 'TET2')), 5)
    expect_equal(length(delete.hypothesis(test_dataset, cause = c('TET2', 'EZH2'))), 5)
    expect_equal(length(delete.hypothesis(test_dataset, effect = c('TET2', 'EZH2'))), 5)
    expect_warning(delete.hypothesis(test_dataset))
    expect_error(delete.hypothesis(test_dataset, event = c('TET2', 'EZH2')))
    expect_error(delete.hypothesis(test_model, event = 'TET2'))
    expect_error(delete.hypothesis(test_dataset, cause = 'XXX'))
    expect_error(delete.hypothesis(test_dataset, effect = 'XXX'))
})

context("delete.pattern")

test_that("delete.pattern is working", {
    expect_equal(length(delete.pattern(test_dataset, pattern = 'XOR_EZH2')), 5)
    expect_error(delete.pattern(test_model, pattern = 'XOR_EZH2'))
    expect_error(delete.pattern(test_dataset, pattern = 'banana'))
    a =  delete.pattern(test_dataset, pattern = 'XOR_EZH2')
    expect_equal(length(delete.pattern(a, pattern = 'OR_CSF3R')), 5)
})

context("delete.model")

test_that("delete.model is working", {
    expect_equal(length(delete.model(test_model)), 5)
    expect_error(delete.model(test_dataset))
})

context("delete.samples")

test_that("delete.model is working", {
    expect_equal(length(delete.samples(test_dataset_no_hypos, samples='patient 1')), 5)
    expect_error(delete.samples(test_dataset_no_hypos, samples='EROS'))
    expect_error(delete.samples(test_model))
})