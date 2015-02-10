# Binds events from one or more datasets, which must be defined over the same set of samples.
# @... the input datasets
ebind = function(...)
{
  # merge two  datasets at a time.
  events.pairwise.bind = function(x, y)
  {  
    is.compliant(x, 'ebind: input x')
    is.compliant(y, 'ebind: input y')
    
    if(!all(as.samples(x) == as.samples(y)))
      stop('Datasets have different samples, won\'t bind!')
    
    z = list()
    
    # Copy genotype matrix, and sets its rownames (samples)
    z$genotypes = cbind(x$genotypes, y$genotypes)
    colnames(z$genotypes) = paste('G', 1:ncol(z$genotypes), sep='')
    
    # Copy annotations for gene symbols etc.
    z$annotations = rbind(x$annotations, y$annotations)
    rownames(z$annotations) = colnames(z$genotypes)
    
    # Copy types
    z$types = unique(rbind(x$types, y$types))	
    
    # Copy stages, if present
    if(has.stages(x) && has.stages(y) && !all(as.stages(x) == has.stages(y), na.rm=TRUE))
      stop('Patients have different stages, won\'t merge!')
    
    if(has.stages(x)) 
    {
      z$stages = as.stages(x)
      colnames(z$stages) = 'stage'
    }
    
    is.compliant(z, 'ebind: output')
    
    return(z)
  }
  
	input = list(...)
	
	if(length(input) <= 1) return(input);
		
	# This could be done with a foldR
	z = events.pairwise.bind(input[[1]], y=input[[2]])
	for(i in 3:length(input))
		z = events.pairwise.bind(z, input[[i]])
		  
	return(z)
 }

# Binds samples from one or more datasets, which must be defined over the same set of events
# @... the input datasets
sbind = function(...)
{
  # merge two datasets at a time.
  samples.pairwise.bind = function(x, y)
  {  
    is.compliant(x, 'sbind: input x')
    is.compliant(y, 'sbind: input y')
    
    if(!all(as.events(x) == as.events(y)))
    {
     # cat(as.events(x)[as.events(x) == as.events(y),])
      stop('Datasets have different events, can not bind!')
    }
    z = list()
    
    # Copy genotypes and annotations
    z$genotypes = rbind(x$genotypes, y$genotypes)
    z$annotations = x$annotations
    z$types = x$types
        
    # Copy stages, if present
    if(has.stages(x) || has.stages(y)) 
    {
      if(has.stages(x)) xstages = as.stages(x)
      else xstages = matrix(rep('NA', nsamples(x)), nrow=nsamples(x))
      
      if(has.stages(y)) ystages = as.stages(y)
      else ystages = matrix(rep('NA', nsamples(y)), nrow=nsamples(y))
      
      z$types = rbind(x$types, y$types)
      
      colnames(z$stages) = 'stage'
    }
    
    
    is.compliant(z, 'sbind: output')
    
    return(z)
  }
  
  input = list(...)
  
  if(length(input) <= 1) return(input);
  
  # This could be done with a foldR
  z = samples.pairwise.bind(input[[1]], y=input[[2]])
  for(i in 3:length(input))
    z = samples.pairwise.bind(z, input[[i]])
  
  return(z)
}