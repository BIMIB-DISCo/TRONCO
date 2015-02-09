# Deletes all events which have frequency 0 in the dataset.
trim = function(x) {
	is.compliant(x, 'trim: input')
  
  x = enforce.numeric(x)
	
	del = names(which(colSums(x$genotypes) == 0))

	x$genotypes = x$genotypes[, !colnames(x$genotypes) %in% del]
	x$annotations = x$annotations[!rownames(x$annotations) %in% del,]

	x$types = matrix(x$types[ unique(x$annotations[,'type']), ncol=1])
	rownames(x$types) = unique(x$annotations[,'type'])
	colnames(x$types) = 'color'
	

	is.compliant(x, 'trim: output')		
	return(x)
}