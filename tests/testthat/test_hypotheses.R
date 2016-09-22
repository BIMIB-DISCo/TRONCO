data(test_dataset_no_hypos)
data(test_model)

context("hypothesis.add")

test_that("hypotheses.add is working", {
    expect_equal(length(hypothesis.add(test_dataset_no_hypos,
        'OR ASXL1 EZH2 TET2',
        OR('ASXL1', 'EZH2', 'TET2'))), 5)
    expect_equal(length(hypothesis.add(test_dataset_no_hypos,
        'XOR ASXL1 EZH2 TET2',
        XOR('ASXL1', 'EZH2', 'TET2'))), 5)
    expect_equal(length(hypothesis.add(test_dataset_no_hypos,
        'AND SETBP1 CSF3R',
        AND('SETBP1', 'CSF3R'))), 5)
    expect_equal(length(hypothesis.add(test_dataset_no_hypos,
        'test',
        OR('ASXL1', 'EZH2', 'TET2'),
        pattern.effect='SETBP1',
        pattern.cause='CSF3R')), 5)
    hypos = hypothesis.add(test_dataset_no_hypos,
        'OR ASXL1 EZH2 TET2',
        OR('ASXL1', 'EZH2', 'TET2'))
    expect_error(hypothesis.add(hypos,
        'OR ASXL1 EZH2 TET2 2',
        OR('ASXL1', 'EZH2', 'TET2')))
    expect_error(hypothesis.add(test_dataset_no_hypos,
        'AND ASXL1 EZH2 TET2',
        AND('ASXL1', 'EZH2', 'TET2')))
    expect_error(hypothesis.add(hypos,
        'test2',
        OR('ASXL1', 'XXX', 'TET2')))
    expect_error(hypothesis.add(hypos,
        'test',
        OR('ASXL1', 'EZH2', 'TET2')))
    expect_error(hypothesis.add(test_model,
        'OR ASXL1 EZH2 TET2 2',
        OR('ASXL1', 'EZH2', 'TET2')))
    expect_error(hypothesis.add(test_dataset_no_hypos,
        'test',
        OR('ASXL1', 'EZH2', 'TET2'),
        pattern.effect='TET2',
        pattern.cause='CSF3R'))
    expect_error(hypothesis.add(test_dataset_no_hypos,
        'test',
        OR('ASXL1', 'EZH2', 'TET2'),
        pattern.effect='XXX',
        pattern.cause='CSF3R'))
    expect_error(hypothesis.add(test_dataset_no_hypos,
        'test',
        OR('ASXL1', 'EZH2', 'TET2'),
        pattern.effect='SETBP1',
        pattern.cause='TET2'))
    expect_error(hypothesis.add(test_dataset_no_hypos,
        'test',
        OR('ASXL1', 'EZH2', 'TET2'),
        pattern.effect='SETBP1',
        pattern.cause='XXX'))
})

context("hypothesis.add.group")

test_that("hypothesis.add.group is working", {
    expect_equal(length(hypothesis.add.group(test_dataset_no_hypos,
        OR,
        c('ASXL1', 'EZH2', 'TET2'),
        silent = TRUE)), 5)
    expect_equal(length(hypothesis.add.group(test_dataset_no_hypos,
        OR,
        c('ASXL1', 'EZH2', 'TET2'),
        min.prob = 0.2,
        silent = TRUE)), 5)
    expect_error(hypothesis.add.group(test_model,
        OR,
        c('ASXL1', 'EZH2', 'TET2'),
        silent = TRUE))
    expect_warning(hypothesis.add.group(test_dataset_no_hypos,
        OR,
        c('ASXL1'),
        silent = TRUE))
    expect_error(hypothesis.add.group(test_model,
        OR,
        c('ASXL1', 'EZH2', 'TET2'),
        dim.max = 2,
        dim.min = 1,
        silent = TRUE))
})

context("hypothesis.add.homologous")

test_that("hypothesis.add.homologous is working", {
    expect_equal(length(hypothesis.add.homologous(test_dataset_no_hypos, silent = TRUE)), 5)
    expect_error(hypothesis.add.homologous(test_model))
})
