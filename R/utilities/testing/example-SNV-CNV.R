# Load Tronco files
source('../loading/import.gistic.R')
source('../loading/import.mutations.R')
source('../loading/import.genotypes.R')
source('../loading/cbio.query.R')
source('../loading/merge.genotypes.R')
source('../loading/union.types.R')
source('../correctness/is.compliant.R')
source('../visualization/oncoprint.R')
source('../visualization/oncoprint.cbio.R')

# Query cbio for GISTIC
genes = c('APC', 'KRAS', 'CTNNB1', 'PTEN', 'PIK3CA', 'ARID1A', 'TP53', 'ATM', 'ERBB2', 'PIK3R1', 'ASXL1', 'BRAF1', 'BRAF2')

cnv.data = cbio.query(22, 11, 9, genes=genes)
mut.data = cbio.query(22, 11, 1, genes=genes)

# Import GISTIC
gistic.bool = import.gistic(cnv.data, merge=F)
mut.bool = import.mutations(mut.data, color='darkgreen')

oncoprint(gistic.bool, ann.stage=F, cellwidth=3, cellheigth=5, null.color='lightgray', hide.zeroes = T, col.cluster=F)
oncoprint(mut.bool, ann.stage=F, cellwidth=3, cellheigth=5, null.color='lightgray', hide.zeroes = T, col.cluster=F)

all = merge.genotypes(mut.bool, gistic.bool)
oncoprint(all, ann.stage=F, cellwidth=3, cellheigth=5, null.color='lightgray', hide.zeroes = T, col.cluster=F)
oncoprint.cbio(all)

# Different ways to see data!
all = union.types(all, 'Homozygous Loss', 'Heterozygous Loss', new.type='Loss', new.color='darkblue')
oncoprint(all, ann.stage=F, cellwidth=3, cellheigth=5, null.color='lightgray', hide.zeroes = T)

all = union.types(all, 'Low-level Gain', 'High-level Gain', new.type='Gain', new.color='darkred')
oncoprint(all, ann.stage=F, cellwidth=3, cellheigth=5, null.color='lightgray', hide.zeroes = T)

all = union.types(all, 'Gain', 'Loss', new.type='CNA', new.color='darkred')
oncoprint(all, ann.stage=F, cellwidth=3, cellheigth=5, null.color='lightgray', hide.zeroes = T)

all = union.types(all, 'CNA', 'SNV', new.type='Alteration', new.color='darkgray')
oncoprint(all, ann.stage=F, cellwidth=3, cellheigth=5, null.color='lightgray', hide.zeroes = T)

