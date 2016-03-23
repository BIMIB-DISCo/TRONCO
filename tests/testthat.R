library(testthat)
library(TRONCO)

data(maf)
muts = import.MAF(maf)
hypo = hypothesis.add(muts, 'test', OR('ABAT', 'ABCC3'))
no_hypo = delete.hypothesis(hypo, 'test')
data(gistic)
gistic = import.GISTIC(gistic)
gistic_model = tronco.caprese(gistic)

test_check("TRONCO")