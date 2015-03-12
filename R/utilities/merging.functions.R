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
    
    if(!(type.one %in% x$annotations[, 'type']) || !(type.two %in% x$annotations[, 'type'])) return(x);
    
    z = list()
    
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
      cat('*** Merge \'', type.one, '\' and \'', type.two, '\'. \n', sep='')
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

      x = rename.type(x, type.one, new.type)
      x = rename.type(x,type.two, new.type)
      x = change.color(x, new.type, new.color)
      z = x
    }
    
    is.compliant(z, 'merge.types: output')
    return(z)
  }
  
  # check if x is compliant
  is.compliant(x)
  
  input = list(...)
  
  if(length(input) == 1 && is.null(input[[1]])) {
    cat(paste('Merging all types: ', paste(as.types(x), collapse=', ', sep=''),' into ', new.type,'.\n', sep=''))
    input = as.list(as.types(x))
  } 

  if(length(input) <= 1) { 
    #cat(paste('Only one type selected: ', input,', done.\n', sep=''))
  
    x = rename.type(x, input, new.type)  
    x = change.color(x, new.type, new.color)
    return(x) 
  }
  
  types.check = lapply(input, function(type){ type %in% as.types(x) })
  if(!all(unlist(types.check))) {
    t = (input[!unlist(types.check)])
    stop('Wrong type \'', t, '\' selected, please check as.types()')
  }
  if(any(duplicated(input))) {
    stop('You tried to merge the same type more than one time, and that\' not a good idea!')
  }
  if(new.type %in% as.types(x)) {
    stop(paste(new.type, 'already in as.types()'))
  }
  
  z = merge.two.types(x, input[[1]], input[[2]], new.type, new.color)
  if (!(length(input) > 2))
    return(z)
  
  for(i in 3:length(input))
    z = merge.two.types(z, new.type, input[[i]], new.type, new.color)
  
  return(z)
  
}



