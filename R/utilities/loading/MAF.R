#
# import.mutations - import from MAF file
#

"import.MAF" <- function(file, sep='\t', is.TCGA=TRUE) 
{
	cat('*** Importing from file: ', file, ' ... ')
  maf = read.delim(file, comment.char = "#", sep = sep, header = TRUE, stringsAsFactors = FALSE)    
  cat('DONE\n')
  
  #### Auxiliary functions to extract information from the MAF file
  # This is the possibly smallest type of information required to prepare a TRONCO file
  # If any necessary information is missing, execution is aborted
  variants = function(x)
  {
    if(!('Variant_Classification' %in% colnames(x)))
      warning('Missing Variant_Classification flag in MAF file.')
    
    return(unique(x[, 'Variant_Classification']))
  }
  
	samples = function(x)
	{
	  if(!('Tumor_Sample_Barcode' %in% colnames(x)))
	    stop('Missing Tumor_Sample_Barcode flag in MAF file - will not import.')
	  
	  return(unique(x[, 'Tumor_Sample_Barcode']))
	}
	
	genes = function(x)
	{
	  if(!('Hugo_Symbol' %in% colnames(x)))
	    stop('Missing Hugo_Symbol flag in MAF file - will not import.')
	  
	  return(unique(x[, 'Hugo_Symbol']))
	}

	valid.calls = function(x)
	{
	  if(!('Validation_Status' %in% colnames(x)))
	    warning('Missing Validation_Status flag in MAF file.')
	  else
	    return(which(x[, 'Validation_Status'] == 'Valid'))
	}
  
  as.TCGA.patients = function(x)
  {
    samples = samples(x)
    patients = substr(samples, 0, 12)  
        
    return(unique(patients))
  }
	
  # General report about the mutations stored in this MAF
	cat('*** MAF report: ')
	if(is.TCGA) cat('TCGA=TRUE')
  
  MAF.variants = variants(maf)
	MAF.samples = samples(maf)
	MAF.genes = genes(maf)
	
  cat('\nType of annotated mutations: \n')
  print(MAF.variants)

  cat('Number of samples:', length(MAF.samples), '\n')

  # If it is TCGA you should check for multiple samples per patient
  if(is.TCGA)
  {
    TCGA.patients = as.TCGA.patients(maf)
    cat('[TCGA = TRUE] Number of TCGA patients:', length(TCGA.patients), '\n')

    if(length(TCGA.patients) != length(MAF.samples))
      warning('This MAF contains duplicate samples for some patients - use TCGA functions for further information')
  }      
  
  which.valid.calls = valid.calls(maf)
  n.valid = length(which.valid.calls)
	cat('Number of annotated mutations:', nrow(maf), '\n')
  
	if(('Validation_Status' %in% colnames(maf)))
	  cat('Mutations annotated with \"Valid\" flag (%):', round(n.valid/nrow(maf)*100, 0), '\n')
	else
	  cat('Mutations annotated with \"Valid\" flag (%): missing flag\n')
	
	cat('Number of genes (Hugo_Symbol):', length(MAF.genes), '\n')
	
	cat('Starting conversion from MAF to 0/1 mutation profiles (1 = mutation) :')
  cat(length(MAF.samples), 'x', length(MAF.genes), '\n')
  
	flush.console()
	pb <- txtProgressBar(1, nrow(maf), style = 3)

  # Temporary binary matrix
  binary.mutations = matrix(0, 
                            nrow=length(MAF.samples), 
                            ncol=length(MAF.genes))
	
  colnames(binary.mutations) = MAF.genes
	rownames(binary.mutations) = MAF.samples
	
	for(i in 1:nrow(maf))
	{  
	  setTxtProgressBar(pb, i)
	  binary.mutations[maf$Tumor_Sample_Barcode[i], maf$Hugo_Symbol[i]] = 1 
	}
	close(pb)
	
	cat('Starting conversion from MAF to TRONCO data type.\n')
	tronco.data = import.genotypes(binary.mutations, event.type = 'Mutation')
  is.compliant(tronco.data)

	return(tronco.data)
}

#data = import.MAF('TCGA_CRC_Suppl_Table2_Mutations_20120719.csv', sep=';', is.TCGA = F)
#show(data)

extract.MAF.HuGO.Entrez.map = function(file, sep='\t') 
{
  cat('*** Importing from file: ', file, ' ... ')
  
  maf = read.delim(file, comment.char = "#", sep = sep, header = TRUE, stringsAsFactors = FALSE)    
  cat('DONE\n')
  
  map = unique(maf[, c('Hugo_Symbol', 'Entrez_Gene_Id')])
  map.missing = which(map[, 'Entrez_Gene_Id'] == 0)
  map.missing = map[map.missing, 'Hugo_Symbol', drop = FALSE]
  
  if(nrow(map.missing) > 0) warning('The are Hugo_Symbol with Entrez_Gene_Id equal to 0')
  else cat('Map seem consistent (non-zero Entrez_Gene_Id) for', nrow(map), 'genes')
  
  return(map)  
}
  