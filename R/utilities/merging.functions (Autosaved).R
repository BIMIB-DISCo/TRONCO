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
	
	print(shared.genes)
	print((tone.genes))
	print(ttwo.genes)
	
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
		z$types = rbind(z$types, c(new.type, new.color))
		rownames(z$types)[length(z$types)] = new.type
	}
	
	# Remove type.one/two from $types flag
	#z$types[type.one, ] = NULL
	#z$types[type.two, ] = NULL


	
	# Now, work with genes which have both types of events
    if(length(shared.genes) > 0) 
    {
    	# These are genotypes restricted to the events we want to process
    	data.tone = as.gene(x, genes = shared.genes, types = type.one)
		data.ttwo = as.gene(x, genes = shared.genes, types = type.two)
				
    	for(i in 1:ncol(data.tone))
      	{
      		# We can delete the genotypes which we just extracted 
      		z = delete.event(z, gene = shared.genes[i], type = type.one)
     		z = delete.event(z, gene = shared.genes[i], type = type.two)
	        
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
	        
	        
	        print(z)
	        print(import.genotypes(geno.tone.ttwo, default.variant = new.type, color = new.color))
	        
	        z = ebind(z, import.genotypes(geno.tone.ttwo, default.variant = new.type, color = new.color))                 
             
                	
        	# Remove the c1/c2 columns from genotypes, and add a new column for c1+c2
        	
            #print(z)		       
            stop('ssssafasafasffss') 
      }
      
      # Restore key access
      colnames(z$genotypes) = paste('G', 1:ncol(z$genotypes), sep='')
	  rownames( z$annotations) = colnames(z$genotypes) 
	}

	#print(z$annotations)
	#print(z$types)

		
		print(z$types) 

		
	stop("z")

	# Change annotations first
	
	# Then check if    
	    
    
    # Extract column names for genotypes
    col.tone = x$annotations[x$annotations[, 'type'] == type.one, ]
    col.tone = cbind(col.tone, refcol=rownames(col.tone))
    col.ttwo = x$annotations[x$annotations[, 'type'] == type.two, ]
    col.ttwo = cbind(col.ttwo, refcol=rownames(col.ttwo))
    
    cols = merge(col.tone, col.ttwo, by='event')
    data.tone = x$genotypes[, as.character(cols[ , 'refcol.x'])]
    data.ttwo = x$genotypes[, as.character(cols[ , 'refcol.y'])]
    z$genotypes = x$genotypes[, !(colnames(x$genotypes) %in% c(colnames(data.tone), colnames(data.ttwo))) ]
    
    un = c(col.tone[, 'event'], col.ttwo[, 'event'])
    
    tone = setdiff( unique(un), cols[,'event'])
    
    z$annotations = x$annotations[colnames(z$genotypes),] 
    
    x = enforce.numeric(x)
    
    # if there is an intersection...
    if (nrow(cols) > 1) {
      # ..merge everything..
      
      # print pb
      	flush.console()

      pb <- txtProgressBar(1, nrow(cols), style = 3)
      
      for(i in 1:nrow(cols))
      {
        c1 = x$genotypes[, as.character(cols[ i, 'refcol.x'])]
        c2 = x$genotypes[, as.character(cols[ i, 'refcol.y'])]
        
        c = c1+c2
        c[c==2] = 1
        
        # Copy genotype matrix and annotations
        z$genotypes = cbind(z$genotypes, c)
        z$annotations = rbind(z$annotations, c(new.type, as.character(cols[ i, 'event'])))
        
        # update pb
        setTxtProgressBar(pb, i)
      }
      
      # close pb
      close(pb)
      
      z$annotations[which(z$annotations[,'type'] == type.one), 'type' ] = new.type
      z$annotations[which(z$annotations[,'type'] == type.two), 'type' ] = new.type
      
      colnames(z$genotypes) = paste('G', 1:ncol(z$genotypes), sep='')
      rownames(z$annotations) = colnames(z$genotypes)
      
      # Copy types
      idx = which(z$annotations[,'type'] != new.type)
      
      z$types = matrix(x$types[ unique(z$annotations[idx,'type'], ) ,])
      rownames(z$types) = unique(z$annotations[idx,'type'], new.type)
      
      z$types = rbind(z$types, new.color)
      rownames(z$types)[length(z$types)] = new.type
      
      colnames(z$types) = 'color'
      
      # Copy stages, if present.  
      if(has.stages(x)) z$stages = x$stages
      
    } else {
      # ..simply rename and change the color
#	 print('Types in: ')	
#	 print(as.types(x))	
#	 print('Types one: ')	
#	 print(type.one)
#	 print('Types new: ')	
#	 print(new.type)
	 	
					
      x = rename.type(x, type.one, new.type)
 #   print('SSS1')
  #  print(x)
    
      x = rename.type(x,type.two, new.type)
   # print('SSS2')

      x = change.color(x, new.type, new.color)
      z = x
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
  {
  	  cat('Merging',  new.type, 'with', input[[i]], 'in', new.type, '\n')

  	  z = merge.two.types(z, new.type, input[[i]], new.type, new.color)
  }
  
  return(z)
  
}



