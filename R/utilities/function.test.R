# load sub1
load('R/utilities/CNV-MUTs-sub1-geno.RData')

# gene test
test = sub1
as.genes(sub1)
test = rename.gene(test, 'APC', 'new name')
as.genes(sub1)
test = delete.gene(test, 'new name')
as.genes(sub1)
is.compliant(test)

# type test
test = sub1
as.types(test)
test = rename.type(test, 'SNV', 'mod')
as.types(test)
test = delete.type(test, 'mod')
as.types(test)
is.compliant(test)

