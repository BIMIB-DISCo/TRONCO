# For an input dataset merge all the events of two distincit types (e.g., say that missense and indel
# mutations are events of a unique "mutation" type)
#
# TODO : extend with ellipsis to allow merging of n types
#
# @ x : the input dataset
# @ type.one: type to merge
# @ type.two: type to merge
# @ new.type: label for the new type to create
# @ new.color: color for the new type to create
merge.types = function(x, type.one, type.two, new.type=paste(type.one, type.two, sep=':'), new.color='khaki')
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
  
	pb = txtProgressBar(1, nrow(cols), style = 3);      
	for(i in 1:nrow(cols))
	{
	  setTxtProgressBar(pb, i)  
	  
		c1 = x$genotypes[, as.character(cols[ i, 'refcol.x'])]
		c2 = x$genotypes[, as.character(cols[ i, 'refcol.y'])]

    c = c1+c2
		c[c==2] = 1
		
		# Copy genotype matrix and annotations
		z$genotypes = cbind(z$genotypes, c)
		z$annotations = rbind(z$annotations, c(new.type, as.character(cols[ i, 'event'])))
	}
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

  is.compliant(z, 'merge.types: output')
	 return(z)
}
