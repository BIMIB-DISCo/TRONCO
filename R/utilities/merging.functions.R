# For an input dataset merge all the events of two or more distincit types (e.g., say that missense and indel
# mutations are events of a unique "mutation" type)
#
# @ x : the input dataset
# @ ...: type to merge
# @ new.type: label for the new type to create
# @ new.color: color for the new type to create
merge.types = function(x, ..., new.type='new.type', new.color='khaki') {
  
  # internal function, merge two types
  merge.two.types = function(x, type.one, type.two, new.type=paste(type.one, type.two, sep=':'), new.color='khaki')
  {  
    is.compliant(x, 'merge.types: input x')
        
    cat('Merging events of type', type.one, 'and', type.two, '-->', new.type, '\n')


    # The returned object
    z = x
    
    # Events of type one/two
    ev.tone = as.events(x, types=type.one)
    ev.ttwo = as.events(x, types=type.two)
    
    # We need the list of genes which have events of both type one and two, and which are specific to a certain type
    shared.genes = intersect(ev.ttwo[, 'event'], ev.tone[, 'event'])
	tone.genes = setdiff(ev.tone[, 'event'], ev.ttwo[, 'event'])
	ttwo.genes = setdiff(ev.ttwo[, 'event'], ev.tone[, 'event'])
	
	# Change names of events which are type specific to type one
	if(length(tone.genes) > 0)
	{
		keys.tone.genes = rownames(as.events(x, genes = tone.genes, types = type.one))
		z$annotations[keys.tone.genes, 'type'] = new.type	   
	}
	
	# Change names of events which are type specific to type two
	if(length(ttwo.genes) > 0)
	{
		keys.ttwo.genes = rownames(as.events(x, genes = ttwo.genes, types = type.two))
		z$annotations[keys.ttwo.genes, 'type'] = new.type	    
	}
	
	# Add new.type to the $types annotation if it does not exist - this leaves "z" still compliant
	if(!(new.type %in% rownames(z$types))) 
	{		
		z$types = rbind(z$types, new.color)
		rownames(z$types)[length(z$types)] = new.type
	}
		
	# Now, work with genes which have both types of events
    if(length(shared.genes) > 0) 
    {
    	# These are genotypes restricted to the events we want to process
    	data.tone = as.gene(x, genes = shared.genes, types = type.one)
		data.ttwo = as.gene(x, genes = shared.genes, types = type.two)

		# We can delete the genotypes which we just extracted 
      	for(i in 1:length(data.tone)) z = delete.event(z, gene = shared.genes[i], type = type.one)
    	for(i in 1:length(data.tone)) z = delete.event(z, gene = shared.genes[i], type = type.two)
     		
		for(i in 1:ncol(data.tone))
      	{
      	    # Get the 1s for each event	        
	        one.tone = data.tone[, i, drop = F] == 1
	        one.ttwo = data.ttwo[, i, drop = F] == 1
	        
	        # Build a simple matrix for that
	        geno.tone.ttwo = matrix(rep(0, nsamples(x)), ncol = 1)
	        geno.tone.ttwo[which(one.tone), ] = 1
	        geno.tone.ttwo[which(one.ttwo), ] = 1
	        
	        # With its colnames 
	        colnames(geno.tone.ttwo) = shared.genes[i]
	        rownames(geno.tone.ttwo) = as.samples(x)    

			# Create TRONCO input and bind datasets	        
	        genotype = import.genotypes(geno.tone.ttwo, event.type = new.type, color = new.color)
	        z = ebind(z, genotype)                
         }
      
  }

	is.compliant(z, 'merge.types (pairwise): output')
    return(z)
  }
  
  # check if x is compliant
  is.compliant(x)
    
  input = list(...)

	# TODO Change this in a better way (deafult ellipsis?)  
  if(length(input) == 1 && is.null(input[[1]])) {
    input = as.list(as.types(x))
  } 

  cat(paste('*** Aggregating events of type(s) {', paste(unlist(input), collapse=', ', sep=''),'} in a unique event with label \"', new.type,'\".\n', sep=''))
  

  if(length(input) <= 1) { 
    cat('One input type provided, using renaming functions.\n')
  
    x = rename.type(x, input, new.type)  
    x = change.color(x, new.type, new.color)
    return(x) 
  }
  
  
  types.check = lapply(input, function(type){ type %in% as.types(x) })
  
  if(!all(unlist(types.check))) {
    t = (input[!unlist(types.check)])
    
    stop('No events of type \'', t, '\' in input dataset, will not merge.')
  }
  
  if(any(duplicated(input))) {
    stop('No duplicated types are allowed, will not merge.')
  }
  
  if(new.type %in% as.types(x)) 
  {
    stop(paste0(new.type, 'is already used in input dataset, will not merge'))
  }
  
  z = merge.two.types(x, input[[1]], input[[2]], new.type, new.color)
  if (!(length(input) > 2))
    return(z)

  
  for(i in 3:length(input))
  	  z = merge.two.types(z, new.type, input[[i]], new.type, new.color)
  
  return(z)
  
}



