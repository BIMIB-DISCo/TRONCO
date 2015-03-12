TCGA.multiple.samples.per.patient = function(x)
{
  is.compliant(x)
  
  samples = as.samples(x)
  samples.truncated = substring(samples, 0, 12)
  
  patients = unique(samples.truncated)
  
  if(length(patients) != nsamples(x))
  {
    dup.samples.start = which(duplicated(samples.truncated)) 
    dup.samples.last = which(duplicated(samples.truncated, fromLast = T))
    
    return(sort(samples[c(dup.samples.start, dup.samples.last)])) 
  }
  else return(NA)
}

TCGA.remove.multiple.samples = function(x)
{
  is.compliant(x, err.fun='Removing TCGA multiple samples (input)')

  dup = TCGA.multiple.samples(x)
  dup.truncated = substring(dup, 0, 12)
  patients = unique(dup.truncated)

  for(i in 1:length(patients))
  {
    patients.samples = which(dup.truncated == patients[i])
    multiple.samples = dup[patients.samples]
    
    cat('Patient', patients[i], 'with sample aliquotes\n' )
    print(substring(multiple.samples, 14, 29))

    keep = max(multiple.samples)
    discard = multiple.samples[which(multiple.samples != keep)]
    
    cat('Selecting', keep, '\n')
    x = delete.samples(x, discard)
  }
    
  is.compliant(x, err.fun='Removing TCGA multiple samples (output)')
  return(x)
}

TCGA.shorten.barcodes = function(x)
{
  is.compliant(x, err.fun='Shartening TCGA barcodes (input)')

  # Check if it has duplicated barcodes
  if(!all(is.na(TCGA.multiple.samples.per.patient(x))))
    stop(
      paste('This dataset contains multiple samples for some patients - cannot consolidate.',
            '\n Samples with barcodes indicating multiple patients: \n', paste(TCGA.multiple.samples.per.patient(x), collapse = '\n'), '.'
            , sep =''))
  
  # Shorten sample barcodes
  rownames(x$genotypes) = substring(rownames(x$genotypes), 0, 12)
  if(has.stages(x)) rownames(x$stages) = rownames(x$genotypes)

  is.compliant(x, err.fun='Shartening TCGA barcodes (output)')
  return(x)    
}
  

TCGA.map.clinical.data = function(file, sep='\t', column.samples, column.map)
{
  
  data = read.delim(
    file = file,
    sep = sep,
    header = TRUE,
    stringsAsFactors=F)
  
  if(!(column.samples %in% colnames(data))) 
    stop(paste('Cannot find samples column \"', column.samples, '\". Available columns: \n\t',
               paste(colnames(data), collapse='\n\t'), sep = ''))
  
  if(!(column.map %in% colnames(data))) 
    stop(paste('Cannot find required map column \"', column.map, '\". Available columns: \n\t',
               paste(colnames(data), collapse='\n\t'), sep = ''))
  
  map = data.frame(data[, column.map], row.names = data[, column.samples])
  colnames(map) = column.map
  
  return(map)
}
