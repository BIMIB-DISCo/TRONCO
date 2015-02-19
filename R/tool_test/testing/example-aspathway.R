# Load Tronco files
source('../visualization/oncoprint.R')
source('../loading/merge.genotypes.R')
source('../pathways/pathways.extract.R')
source('../selection/events.selection.R')

# Use only cluster 1, available in tcga.sub1
source('./example-dataloading.R')

# Pathway members according to TCGA paper
Wnt = c('APC', 'CTNNB1', 'DKK1', 'DKK2', 'DKK3', 'DKK4', 'LRP5', 'FZD10', 'FAM123B', 'AXIN2', 'TCF7L2', 'FBXW7', 'ARID1A', 'SOX9')
RAS = c('ERBB2', 'ERBB3', 'NRAS', 'KRAS', 'BRAF')
PI3K = c('IGF2', 'IRS2', 'PIK3CA', 'PIK3R1', 'PTEN')
TGFb = c('TGFBR1', 'TGFBR2', 'ACVR1B', 'ACVR2A', 'SMAD2', 'SMAD3', 'SMAD4' )
P53 = c('TP53', 'ATM')   

# Extract a pathway and visualize it (annotations disabled)
p.pi3k = as.pathway(tcga.sub1, PI3K, 'PI3K')
oncoprint(p.pi3k)

# Extract othe pathways
p.wnt = as.pathway(tcga.sub1, Wnt, 'Wnt')
p.ras = as.pathway(tcga.sub1, RAS, 'RAS')
p.tgfb = as.pathway(tcga.sub1, TGFb, 'TGFb')
p.p53 = as.pathway(tcga.sub1, P53, 'P53')

# # Merge  pathways and visualize them (annotations disabled)
z = merge.genotypes(p.pi3k, p.wnt, p.ras, p.tgfb, p.p53)
oncoprint(z) 



# p.wnt = as.pathway(tcga.sub2, Wnt, 'Wnt')
# p.p53 = as.pathway(tcga.sub2, P53, 'P53')
# z = merge.genotypes(p.wnt, p.p53)
# oncoprint(z) 

