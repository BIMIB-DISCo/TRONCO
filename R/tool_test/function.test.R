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

# ebind test
t1 = sub1
t1 = delete.type(t1, 'Amplification')
t1 = delete.type(t1, 'SNV')
t2 = sub1
t2 = delete.type(t2, 'Deletion')
test = ebind(t1, t2)
is.compliant(test)

# sbind test
t1 = sub1
t1$genotypes = t1$genotypes[1:20,]
t1$stages = t1$stages[ which(rownames(t1$stages) %in% rownames(t1$genotypes)), , drop=F]
is.compliant(t1)
t2 = sub1
t2$genotypes = t2$genotypes[21:40,]
t2$stages = t2$stages[ which(rownames(t2$stages) %in% rownames(t2$genotypes)), , drop=F]
is.compliant(t2)
test = sbind(t1, t2)
is.compliant(test)
