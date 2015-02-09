trim = function(x) {
	is.compliant(x, 'trim: input')
	
	del = names(which(colSums(x$genotypes) == 0))

	x$genotypes = x$genotypes[, !colnames(x$genotypes) %in% del]
	x$annotations = x$annotations[!rownames(x$annotations) %in% del,]

	# print(!colnames(x$genotypes) %in% del)
	# print(!rownames(x$annotations) %in% del)
	# print(x$types)
	# # print(unique(x$annotations[,'type']))
	 # print(str(x$types[ unique(x$annotations[,'type']), ]))
	# print(str(x))
	
	# print(class(x$types))
	
	
	
	x$types = matrix(x$types[ unique(x$annotations[,'type']), ncol=1])
	rownames(x$types) = unique(x$annotations[,'type'])
	colnames(x$types) = 'color'

	
# #   	print(x$types)
	# print(t(x$types))
	
	
	# print(str(x))
	# print(class(x$types))

	is.compliant(x, 'trim: output')		
	return(x)
}