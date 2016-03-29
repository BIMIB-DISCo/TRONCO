data(test_dataset_no_hypos)

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
})

context("hypothesis.add.group")

test_that("hypothesis.add.group is working", {
    expect_equal(length(hypothesis.add.group(test_dataset_no_hypos,
        OR,
        c('ASXL1', 'EZH2', 'TET2'),
        silent = TRUE)), 5)
})

context("hypothesis.add.homologous")

test_that("hypothesis.add.homologous is working", {
    expect_equal(length(hypothesis.add.homologous(test_dataset_no_hypos, silent = TRUE)), 5)
})