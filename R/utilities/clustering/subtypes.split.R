subtypes.split <- function(x, clusters, idx=NA) 
{
	is.compliant(x)
	
	data = x$genotypes

	cat('*** Splitting genotypes according to clustering indexes\n')

	# TODO: check that clusters has at least 2 columns
	# if(!is.matrix(data) != nrow(clusters)) 
		# stop("Error: no concordance among number of samples and clustering assignement.");

	if(nrow(data) != nrow(clusters)) 
		stop("Error: no concordance among number of samples and clustering assignement.");

	if(!all(rownames(clusters) %in% rownames(data)))
		stop('Not all samples are assigned to a cluster!')

	cluster.labels = unique(clusters)
	num.clusters = nrow(cluster.labels)
	
	if(!is.na(idx))
	{
		y = list()
		samples.in.cluster = rownames(clusters)[clusters == idx]
		
		cat(paste('Cluster labeled ', idx, ' has ', length(samples.in.cluster), ' components.\n', sep=''))
		cat(paste('Returning only cluster with label idx=', idx, '.\n', sep=''))
		
		y$genotypes = data[samples.in.cluster, ];
		y$annotations = x$annotations
		y$types = x$types
				
		if(!is.null(x$stages)) 
		{
			y$stages = as.matrix(x$stages[samples.in.cluster, ]) 
			rownames(y$stages) = samples.in.cluster
			# print(head(y$stages))
			# print(head(x$stages))

			# y$stages = x$stages[samples.in.cluster, ]
		}
		
		is.compliant(y, 'subtypes.split with index', stage=!is.null(y$stages))
		return(y)	
		
		
	}
	
	partitions = list()
	for (i in 1:num.clusters) 
	{
		y = list()
		samples.in.cluster = rownames(clusters)[clusters == cluster.labels[i,1]]
		
		cat(paste('Cluster labeled ', cluster.labels[i,1], ' has ', length(samples.in.cluster), ' components:', sep=''))
		
		y$genotypes = data[samples.in.cluster, ];
		y$annotations = x$annotations
		y$types = x$types
				
		if(!is.null(x$stages)) {
			y$stages = as.matrix(x$stages[samples.in.cluster, ]) 
			rownames(y$stages) = samples.in.cluster
			}
		
		is.compliant(y, 'subtypes.split partitionig', stage=!is.null(y$stages))
		partitions = c(partitions, y)	

		cat(' OK.\n')
	}
		
	return(partitions)
}