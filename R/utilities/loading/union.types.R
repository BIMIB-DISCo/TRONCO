#
# union.types - merge two input genotypes if they have same sample set, and stages if any. No duplicates are removed.
#
union.types = function(x, type.one, type.two, new.type=paste(type.one, type.two, sep=':'), new.color='khaki')
{	
	is.compliant(x, 'union.types: input x')
	
	if(!(type.one %in% x$annotations[, 'type']) || !(type.two %in% x$annotations[, 'type'])) return(x);
	
	# TODO check that events are the same and mergeable
	
	z = list()

	# Extract column names for genotypes
	col.tone = x$annotations[x$annotations[, 'type'] == type.one, ]
	col.tone = cbind(col.tone, refcol=rownames(col.tone))
	col.ttwo = x$annotations[x$annotations[, 'type'] == type.two, ]
	col.ttwo = cbind(col.ttwo, refcol=rownames(col.ttwo))
	
	# cat(paste('#', ncol())
	
	cols = merge(col.tone, col.ttwo, by='event')
	data.tone = x$genotypes[, as.character(cols[ , 'refcol.x'])]
	data.ttwo = x$genotypes[, as.character(cols[ , 'refcol.y'])]
	z$genotypes = x$genotypes[, !(colnames(x$genotypes) %in% c(colnames(data.tone), colnames(data.ttwo))) ]

	un = c(col.tone[, 'event'], col.ttwo[, 'event'])
	
	# print(nrow(col.tone))
	# print(nrow(col.ttwo))
	# print(un)
	tone = setdiff( unique(un), cols[,'event'])
	
	
	# print(length(tone))
	# print(unique (c(col.tone[1:10, 'event'], col.ttwo[1:8, 'event'])))
	# print(col.tone)	
	
	z$annotations = x$annotations[colnames(z$genotypes),] 

	for(i in 1:nrow(cols))
	{
		c1 = x$genotypes[, as.character(cols[ i, 'refcol.x'])]
		c2 = x$genotypes[, as.character(cols[ i, 'refcol.y'])]

		# c = c1
		# for(j in 1:length(c1))
			# c[j] = (c1[j] == 1 || c2[j] == 1)
		c = c1+c2
		c[c==2] = 1
		
		# Copy genotype matrix and annotations
		z$genotypes = cbind(z$genotypes, c)
		z$annotations = rbind(z$annotations, c(new.type, as.character(cols[ i, 'event'])))
	}
	
		# print(
	 # which(z$annotations[,'type'] == type.one )
	# )
	
	
	
	# print(colnames(z$genotypes))

	z$annotations[which(z$annotations[,'type'] == type.one), 'type' ] = new.type
	z$annotations[which(z$annotations[,'type'] == type.two), 'type' ] = new.type


	colnames(z$genotypes) = paste('G', 1:ncol(z$genotypes), sep='')
	rownames(z$annotations) = colnames(z$genotypes)

	# print(z$annotations)

	# Copy types
	idx = which(z$annotations[,'type'] != new.type)
	
	# print(unique(z$annotations[idx,'type'], ))
	
	
	z$types = matrix(x$types[ unique(z$annotations[idx,'type'], ) ,])
	rownames(z$types) = unique(z$annotations[idx,'type'], new.type)
	
	z$types = rbind(z$types, new.color)
	rownames(z$types)[length(z$types)] = new.type
	
	colnames(z$types) = 'color'
	
	
	# print(rownames(x$types))
	# print(rownames(z$types))	
	
	# print(z)
	
 	# Copy stages, if present. test compliance	
	if(!is.null(x$stages)) 
	{
		z$stages = x$stages
		is.compliant(z, 'union.types: output', TRUE)
	}
	else is.compliant(z, 'union.types: output')


	 return(z)
}
