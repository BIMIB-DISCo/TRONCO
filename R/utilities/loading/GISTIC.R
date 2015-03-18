#
# import.gistic - Convert a GISTIC score file to TRONCO input.
#


import.GISTIC <- function(x, stage.annot = NA) {
  
  cat('*** GISTIC input format conversion started.\n'); 	
  
  # For next operations it is convenient to have everything as 'char' rather than 'int'
  if(typeof(x[,1]) != typeof("somechar"))
  {
    cat('Converting input data to character for import speedup.\n')
    rn = rownames(x)
    x = apply(x, 2, as.character) 
    rownames(x) = rn		
  }	
  
  if(any(is.na(x))) warning('NA entries were replaced with 0s.\n')
  x[is.na(x)] = 0		
  
  cat('Creating ', 4 * (ncol(x)), ' events for ', ncol(x) ,'genes \n');
  
  # gene symbols
  enames <- colnames(x)
  if(is.null(enames)) stop('Error: gistic file has no column names and can not imported, aborting!')
  
  cat('Extracting \"Homozygous Loss\" events (GISTIC = -2) \n'); 	
  d.homo = x
  d.homo[d.homo != -2] <- 0
  d.homo[d.homo == -2] <- 1
  
  cat('Extracting \"Heterozygous Loss\" events (GISTIC = -1) \n');   
  d.het  <- x
  d.het[d.het != -1] <- 0
  d.het[d.het == -1] <- 1
  
  cat('Extracting \"Low-level Gain\" events (GISTIC = +1) \n');   
  d.low  <- x
  d.low[d.low != 1] <- 0
  
  cat('Extracting \"High-level Gain\" events (GISTIC = +2) \n');   
  d.high <- x
  d.high[d.high != 2] <- 0
  d.high[d.high == 2] <- 1
    
  cat('Transforming events in TRONCO data types ..... \n')
  d.homo = trim(import.genotypes(d.homo, event.type='Homozygous Loss', color = 'dodgerblue4'))
  d.het  = trim(import.genotypes(d.het, event.type='Heterozygous Loss', color = 'dodgerblue1'))
  d.low  = trim(import.genotypes(d.low, event.type='Low-level Gain', color = 'firebrick1'))
  d.high  = trim(import.genotypes(d.high, event.type='High-level Gain', color = 'firebrick4'))
  
  d.cnv.all = ebind(d.homo, d.het, d.low, d.high)
  
  cat('*** Data extracted, returning only events observed in at least one sample \n', 
      'n=', nevents(d.cnv.all), '(events)\n',
      '|G|=', ngenes(d.cnv.all),'(genes)\n', 
      'm=', nsamples(d.cnv.all), ' (samples)\n')
  
  
  is.compliant(d.cnv.all, 'import.gistic: output')
  return(d.cnv.all);
}
