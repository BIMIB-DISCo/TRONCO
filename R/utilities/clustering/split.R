# TODO: check
# Split cohort (samples) into groups, return either all groups or a specific group.
#
# x: cohort
# clusters: groups, a map  for each sample in 'x' to a group.  
ssplit <- function(x, clusters, idx=NA) 
{
	is.compliant(x)
	
	data = x$genotypes

	cat('*** Splitting cohort into groups.\n')

	# TODO: check that clusters has at least 2 columns
	# if(!is.matrix(data) != nrow(clusters)) 
		# stop("Error: no concordance among number of samples and clustering assignement.");

  # Check that map has correct size
	if(nsamples(x) != nrow(clusters)) 
		stop(paste("Error: cannot split, number of samples (", nsamples(x) , 
               ") and groups (", nrow(clusters),") do not match.", sep=''));

	# Check it is an actual map
	if(!all(rownames(clusters) %in% rownames(data)))
    stop(paste('Error: samples', paste(
      rownames(clusters)[!rownames(clusters) %in% as.samples(x)],
      collapse=', ', sep='')  ,'are not assigned to a group', sep=''));

  # Groups info
  cluster.labels = unique(clusters)
	num.clusters = nrow(cluster.labels)
	
  # Extract a specific group
	if(!is.na(idx))
	{
		y = list()
		samples.in.cluster = rownames(clusters)[clusters == idx]
		
		cat(paste('Group \"', idx, '\" has ', length(samples.in.cluster), 
		          ' samples, returning this group.\n', sep=''))
		
		y$genotypes = data[samples.in.cluster, ];
		y$annotations = x$annotations
		y$types = x$types
				
		if(!is.null(x$stages)) 
		{
			y$stages = as.matrix(x$stages[samples.in.cluster, ]) 
			rownames(y$stages) = samples.in.cluster
		}
		
		is.compliant(y, 'ssplit.split with index')
		return(y)		
	}
	
	# Extract all groups
	partitions = list()
	for (i in 1:num.clusters) 
	{
		y = list()
		samples.in.cluster = rownames(clusters)[clusters == cluster.labels[i,1]]
		
		cat(paste('Group \"', idx, '\" has ', length(samples.in.cluster), 
		          ' samples.\n', sep=''))
		
		y$genotypes = data[samples.in.cluster, ];
		y$annotations = x$annotations
		y$types = x$types
				
		if(!is.null(x$stages)) {
			y$stages = as.matrix(x$stages[samples.in.cluster, ]) 
			rownames(y$stages) = samples.in.cluster
			}
		
		is.compliant(y, 'subtypes.split partitionig')
		partitions = c(partitions, y)	
	}
		
	return(partitions)
}

# Split events into groups according to their types.
#
# x: cohort
tsplit <- function(x) 
{
# Parametro per estrarre un tipo solo di evento?
}