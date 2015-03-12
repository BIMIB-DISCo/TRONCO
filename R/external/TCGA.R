TCGA.multiple.samples = function(x)
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

  dup = TCGA.duplicated.samples(x)
  dup.truncated = substring(dup, 0, 12)
  patients = unique(dup.truncated)

  for(i in 1:length(patients))
  {
    cat('Patient', patients[i], 'with sample aliquotes\n' )
    patients.samples = which(dup.truncated == patients[i])
    print(substring(dup[patients.samples], 14, 29))
    cat('Selecting', max(dup[patients.samples]), '\n')
  }
    
  is.compliant(x, err.fun='Removing TCGA multiple samples (output)')
  return(x)
}

TCGA.shorten.barcodes = function(x)
{
  is.compliant(x, err.fun='Shartening TCGA barcodes (input)')

  # Check if it has duplicated barcodes
  if(!all(is.na(TCGA.duplicated.samples(x))))
    stop(
      paste('This dataset contains multiple samples for some patients - cannot consolidate.',
            '\n Samples with barcodes indicating multiple patients: \n', paste(TCGA.duplicated.samples(x), collapse = '\n'), '.'
            , sep =''))
  
  # Shorten sample barcodes
  rownames(x$genotypes) = substring(rownames(x$genotypes), 0, 12)
  if(has.stages(x)) rownames(x$stages) = rownames(x$genotypes)

  is.compliant(x, err.fun='Shartening TCGA barcodes (output)')
  return(x)    
}
  
  