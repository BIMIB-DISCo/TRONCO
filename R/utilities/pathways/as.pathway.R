# TODO refactor required
#
# as.pathway - Given a cohort and a pathway, return the cohort with events restricted to genes 
# involved in the pathway. This might contain a new 'pathway' genotype with an alteration mark if
# any of the involved genes are altered.
#
# @pathway.genes: gene - symbols - involved in the pathway
# @pathway.name: pathway name
# @pathway.color: pathway color (for visualization)
# @aggregate.pathway: if TRUE does not show the genes in the pathway
as.pathway <- function(x, pathway.genes, pathway.name, 
                       pathway.color='darkgreen', aggregate.pathway = TRUE) 
{
	is.compliant(x, 'as.pathway: input')
	
	data = x$genotypes

	cat(paste('*** Extracting events for pathway: ', pathway.name,'.\n', sep=''))
	
	# Select only those events involving a gene in pathway.genes which is also in x
	y = events.selection(x, NA, filter.in.names=pathway.genes, NA)
	
	# Extend genotypes
  y = enforce.numeric(y)
	
	pathway.genotype = rowSums(y$genotypes)
	pathway.genotype[pathway.genotype > 1] = 1 # Any hit is enough

  # This implemetation is ugly - TODO refactor as follows:
  # - z = import.genotypes(pathway.genotype, ...) # import 
  # - then if agg=F use gbind(y, z) else return z.
  
	y$genotypes = cbind(y$genotypes, pathway.genotype)
	colnames(y$genotypes)[ncol(y$genotypes)] = pathway.name
	y$annotations = rbind(y$annotations, c('Pathway', pathway.name))
	rownames(y$annotations)[nrow(y$annotations)] = pathway.name	
	
	y$types = unique(rbind(y$types, pathway.color))
	rownames(y$types)[nrow(y$types)] = 'Pathway'

	if(aggregate.pathway)
	{
		idx = which(y$annotations[, 'type'] != 'Pathway', arr.ind=T)
		idx.in = which(y$annotations[, 'type'] == 'Pathway', arr.ind=T)

		y$annotations = matrix(y$annotations[idx.in, ], ncol=2)
		rownames(y$annotations) = pathway.name

		y$types = data.frame(y$types['Pathway', ], row.names = 'Pathway')
		# colnames(y$types) = 'color'
		
		# print(y$annotations)
		# print(y$genotypes[, colnames(y$genotypes) %in% names(idx.in)])
		y$genotypes = data.frame(y$genotypes[, colnames(y$genotypes) %in% names(idx.in)])
		colnames(y$genotypes) = pathway.name

		# y = import.genotypes( matrix(pathway.genotype, ncol=1), default.variant='Pathway', color='khaki')
		 # y$genotypes = pathway.genotype
		# colnames(y$genotypes)[ncol(y$genotypes)] = pathway.name
		# y$annotations = rbind(y$annotations, c('pathway', pathway.name))
		# rownames(y$annotations)[nrow(y$annotations)] = pathway.name
		# y$types = unique(rbind(y$types, c(pathway.color)))
		# rownames(y$types)[nrow(y$types)] = 'pathway'
	}
	
	colnames(y$types) = c('color')
	colnames(y$annotations) = c('type', 'event')

	if(!is.null(x$stages))
	{	
		
		y$stages = x$stages
		colnames(y$stages) = 'stage'
	}
# 				print(y)

	is.compliant(y, 'as.pathway: output')
	

	return(y)
}