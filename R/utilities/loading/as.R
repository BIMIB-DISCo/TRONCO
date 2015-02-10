
# Return all samples IDs in the cohort
as.samples = function(x)
{
	return(rownames(x$genotypes))
}

# Return alla gene symbols in the cohort
as.genes = function(x)
{
	return(unlist(unique(x$annotations[, 'event'])))
}

# Return all events for each gene in the cohort
as.events = function(x)
{
  return(x$annotations[,c('type', 'event')])
}

# Return alla stages per sample, if present
as.stages = function(x)
{
  if(has.stages(x))
    return(x$stages)
  else return(NA)
}

# Return all types (labels) of events in the dataset
as.types = function(x)
{
  return(rownames(x$types))
}

# Return all colors
as.colors = function(x)
{
  return(x$types[, 'color'])
}


# Return genotypes for a specific gene (uses all events)
as.gene = function(x, g)
{
	data = data.frame(x$genotypes[, which(as.events(x)[, 'event'] == g)], row.names = as.samples(x))
	colnames(data) = x$annotations[which(as.events(x)[, 'event'] == g)]
	
	return(data)
}

# Return true if stages are present
has.stages = function(x)
{
  return(!is.null(x$stages))
}


# Return true if duplicate events are present
has.duplicates = function(x) {
  # find duplicate over the dataset
  dup = duplicated(as.events(x))
  
  # return true if al least one duplicate is found
  return(any(dup))
}

duplicates = function(x) {
  as.events(x)[duplicated(as.events(x)),]
}


# Short report for a dataset 
show = function(x, view = 10)
{
	is.compliant(x)
	
	cat(paste('Dataset: n=', nsamples(x), ', m=', nevents(x), ', |G|=', ngenes(x), '.\n', sep=''))
	cat(paste('Events: ', paste(as.types(x), collapse=', '), '.\n', sep=''))
	cat(paste('Colors: ', paste(as.colors(x), collapse=', '), '.\n', sep=''))

	if(has.stages(x))
	{
		cat(paste('Stages: '))
s		
		s = unlist(unique(as.stages(x)[, 1]))
		cat(paste(s, collpase=', '), '.\n', sep='')
	 }

	cat(paste('Events: ', view, ' shown.\n'))
	cat(paste(paste('\t', rownames(as.events(x)[1: view,]), ':', as.events(x)[1: view, 1], as.events(x)[1: view, 2], sep=' '), collapse='\n'))

	cat(paste('\nGenotypes for ', view, ' events.\n'))
	head(x$genotypes[,1:view])
}


nsamples = function(x)
{
	return(nrow(x$genotypes))
}

nevents = function(x)
{
	return(ncol(x$genotypes))
}

ngenes = function(x)
{
	return(length(as.genes(x)))
}

# Check internal types and convert accordingly
enforce.numeric = function(x)
{
  if(is.numeric(x$genotypes[1,1]))
  {
    rn = as.samples(x)
    x$genotypes = apply(x$genotypes, 2, as.numeric)
    rownames(x$genotypes) = rn
  }
  
  return(x)
}

enforce.string = function(x)
{
  if(is.character(x$genotypes[1,1]))
  {
    rn = as.samples(x)
    x$genotypes = apply(x$genotypes, 2, as.character)
    rownames(x$genotypes) = rn
  }

  return(x)
}

