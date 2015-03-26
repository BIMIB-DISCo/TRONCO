# Binds events from one or more datasets, which must be defined over the same set of samples.
# @... the input datasets
ebind = function(...)
{
  # merge two  datasets at a time.
  events.pairwise.bind = function(x, y)
  {  
    is.compliant(x, 'ebind: input x')
    is.compliant(y, 'ebind: input y')
	
    samples.intersect = intersect(as.samples(x), as.samples(y))
    if(!(setequal(samples.intersect, as.samples(x)) && setequal(samples.intersect, as.samples(y))))
      stop('Datasets have different samples, won\'t bind!')
    
    z = list()

	y$genotypes = y$genotypes[rownames(x$genotypes), , drop = FALSE]
    y$stages = y$stages[rownames(y$genotypes), , drop=FALSE]
     
    # Copy genotype matrix, and sets its rownames (samples)
    z$genotypes = cbind(x$genotypes, y$genotypes)
    colnames(z$genotypes) = paste('G', 1:ncol(z$genotypes), sep='')
    
    # Copy annotations for gene symbols etc.
    z$annotations = rbind(x$annotations, y$annotations)
    rownames(z$annotations) = colnames(z$genotypes)
    
    # Copy types
    z$types = unique(rbind(x$types, y$types))	

    # print(x)
    # print(y)
    # print(has.stages(x))
    # print(has.stages(y))
    # print(as.stages(x) == as.stages(y))
    
    # Copy stages, if present 
    if(has.stages(x) && has.stages(y))
    {
    		stages.x = as.stages(x)
    		stages.y = as.stages(y)
    		
    		stages.x = stages.x[!is.na(stages.x)]
    	stages.y = stages.y[!is.na(stages.y)]
    	
     	if(any(stages.x != stages.y))
      		stop('Patients have different stages, won\'t merge!')
    }
    
    if(has.stages(x)) 
    		z = annotate.stages(z, as.stages(x))

    is.compliant(z, 'ebind: output')
    
    return(z)
  }
  
  
  input = list(...)

  cat('*** Binding events for', length(input), 'datasets.\n')
  return(Reduce(events.pairwise.bind, input))
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
      
      z$stages = (rbind(x$stages, y$stages))
      
      colnames(z$stages) = 'stage'
    }
    
    
    is.compliant(z, 'sbind: output')
    
    return(z)
  }
  
  input = list(...)
  return(Reduce(samples.pairwise.bind, input))
}