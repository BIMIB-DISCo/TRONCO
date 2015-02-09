#
# aux.pairwise.merge - merge two input genotypes if they have same sample set, and stages if any. No duplicates are removed.
#
aux.pairwise.merge = function(x, y)
{	
	is.compliant(x, 'merge.genotypes: input x')
	is.compliant(y, 'merge.genotypes: input y')
	
	if(!all(rownames(x$genotypes) == rownames(y$genotypes)))
		stop('Genotypes have different samples, won\'t merge!')
	
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
	# print(is.null(x$stages))
	# print(is.null(y$stages))
	# print(head(x$stages))
	# print(head(y$stages))
	# print(x$stages == y$stages)
	if(!(is.null(x$stages)) && !is.null(y$stages) &&
		!all(x$stages == y$stages, na.rm=TRUE))
		stop('Patients have different stages, won\'t merge!')
	
	if(!is.null(x$stages)) z$stages = x$stages
	if(!is.null(y$stages)) z$stages = y$stages
	if(!is.null(z$stages)) 
	{
		colnames(z$stages) = 'stage'
		is.compliant(z, 'merge.genotypes: output', TRUE)
	}
	else is.compliant(z, 'merge.genotypes: output')

	return(z)
}

#
# merge.genotypes - merge a list of input genotypes if they have same sample set, and stages if any. No duplicates are removed.
#
merge.genotypes = function(...)
{
	input = list(...)
	
	if(length(input) <= 1) return(NULL);
		
	# This could be done with a foldR
	z = aux.pairwise.merge(input[[1]], y=input[[2]])
	if(length(input) == 2) return(z)
		
	for(i in 3:length(input))
		z = aux.pairwise.merge(z, input[[i]])
		
	return(z)
 }