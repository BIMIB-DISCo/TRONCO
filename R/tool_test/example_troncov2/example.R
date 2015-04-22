############################################################ LOAD THE LIBRARY
library("TRONCO")
############################################################


###################################################### Somatic Mutations
# TCGA MAF file - originally Excel file, now exported in csv format (separator ";")
# download: 12 March 2015
# URL: https://tcga-data.nci.nih.gov/docs/publications/coadread_2012/

MAF.file = './TCGA_CRC_Suppl_Table2_Mutations_20120719.csv'

# Transform the MAF dataset in a TRONCO (compliant) object
MAF = import.MAF(file = MAF.file, is.TCGA = TRUE, sep = ';')
MAF = change.color(MAF, 'Mutation', 'green')
show(MAF)

###################################################### Check for duplicated samples in the MAF
# Check for duplicated samples - we find them
TCGA.multiple.samples(MAF)

# Remove duplicated samples according to TCGA criteria, add stages - the HuGO map does not change
MAF.patients = TCGA.remove.multiple.samples(MAF)
MAF.patients = TCGA.shorten.barcodes(MAF.patients)

# Save the MAF as a TRONCO compliant Rdata
save(MAF.patients, file = './MAF.patients.Rdata')
show(MAF.patients)

# Selection of mutations with hit-rate > 15%
# Added CTNNB1 
# Removed TTN
MAF.selection = events.selection(MAF.patients, filter.freq=.15, filter.in.names=c('CTNNB1', 'NRAS'), filter.out.names='TTN')
oncoprint(MAF.selection)

# Added hypotheses APC or CTNNB1 - KRAS xor NRAS
data = hypothesis.add(MAF.selection, 'APC or CTNNB1', OR('APC', 'CTNNB1'), '*')
data = hypothesis.add(data, 'KRAS xor NRAS', XOR('KRAS', 'NRAS'), '*')

# Reconstruct with the defaults (BIC as regularization)
crc.bic = tronco.capri(data)
tronco.plot(crc.bic, hidden.and=F)

# Reconstruct with AIC as regularization
crc.aic = tronco.capri(data, REGULARIZATION='aic')
tronco.plot(crc.aic, hidden.and=F)

# Plot both the previous two reconstruction in one unique plot
tronco.plot(crc.bic, secondary=crc.aic, hidden.and=F)
