# 
# as.pathway - Given a genotype 'x' and a pathway, return the subset of 'x' restricted to genes 
# involved in the pathway, plus a new genotype for the pathway which contains an alteration if
# any of the involved genes are altered.
#
# @pathway.genes: gene symbols involved in the pathway
# @pathway.name: pathway symbol
# @pathway.color: pathway color (for visualization)
as.pathway <- function(x, pathway.genes, pathway.name, pathway.color='darkgreen') 
{
	is.compliant(x, 'as.pathway: input')
	
	data = x$genotypes

	cat(paste('*** Extracting events for pathway: ', pathway.name,'.\n', sep=''))
	
	# Select only those events involving a gene in pathway.genes which is also in x
	
	# Select data
	y = events.selection(x, NA, filter.in.names=pathway.genes, NA)
	
	# Extend genotypes
	pathway.genotype = rowSums(y$genotypes)
	pathway.genotype[pathway.genotype>1] = 1 # Any hit is enough
	
	y$genotypes = cbind(y$genotypes, pathway.genotype)
	colnames(y$genotypes)[ncol(y$genotypes)] = pathway.name

	y$annotations = rbind(y$annotations, c('pathway', pathway.name))
	rownames(y$annotations)[nrow(y$annotations)] = pathway.name
	
	y$types = unique(rbind(y$types, c(pathway.color)))
	rownames(y$types)[nrow(y$types)] = 'pathway'
	
	if(!is.null(x$stages))
	{	
		y$stages = x$stages
		colnames(y$stages) = 'stage'
	}
	
	is.compliant(y, 'as.pathway: output')
		
	return(y)
}