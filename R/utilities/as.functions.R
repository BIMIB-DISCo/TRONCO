
# Return all samples IDs in the cohort.
#
# @x: the dataset.
# @genes: a list of genes to query, if NA all genes are used.
as.samples = function(x)
{
	return(rownames(x$genotypes))
}

# Return all gene symbols for which a certain type of event exists.
#
# @x: the dataset.
# @types: the types of events to consider, if NA all available types are used.
as.genes = function(x, types=NA)
{
	return(unlist(unique(as.events(x, types=types)[, 'event'])))
}

# Return all events involving certain genes and types.
#
# @x: the dataset.
# @genes: a list of genes to query, if NA all genes are used.
# @types: the types of events to consider, if NA all available types are used.
as.events = function(x, genes=NA, types=NA)
{
  ann = x$annotations[, c('type', 'event')]

  if(!any(is.na(genes))) ann = ann[ which(ann[, 'event'] %in% genes) , , drop=FALSE] 
  if(!any(is.na(types))) ann = ann[ which(ann[, 'type'] %in% types) , , drop=FALSE]   

  return(ann)
}

# Return all stages per sample, if present.
#
# @x: the dataset.
as.stages = function(x)
{
  if(has.stages(x)) return(x$stages)
  else return(NA)
}

# Return the (labels) of events in the dataset for a set of genes.
#
# @x: the dataset.
# @genes: a list of genes to query, if NA all available genes are used.
as.types = function(x, genes=NA)
{
  return(unlist(unique(as.events(x, genes=genes)[, 'type'])))
}

# Return all colors map for the events in the dataset.
#
# @x: the dataset.
as.colors = function(x)
{
  return(x$types[, 'color'])
}

# Return genotypes for a set of genes and events. 
#
# @x: the dataset.
# @genes: a list of genes to query.
# @types: the types of events to consider, if NA all available types are used.
as.gene = function(x, genes, types=NA)
{
  keys = as.events(x, genes=genes, types=types)
  
	data = data.frame(x$genotypes[, rownames(keys)], row.names = as.samples(x))
  colnames(data) = keys[, 'type']
	
	return(data)
}

# Return true if stages are present
#
# @x: the dataset.
has.stages = function(x)
{
  return(!is.null(x$stages))
}


# Return true if duplicate events are present.
#
# @x: the dataset.
has.duplicates = function(x) {

  # find duplicate over the dataset
  dup = duplicated(as.events(x))
  
  # return true if at least one duplicate is found
  return(any(dup))
}

# Return duplicated events, if any.
#
# @x: the dataset.
duplicates = function(x) {
  as.events(x)[duplicated(as.events(x)),]
}


# Print a short report of a dataset. 
#
# @x: the dataset.
# @view: the number of events to show.
show = function(x, view = 10)
{
	is.compliant(x)
  x = enforce.numeric(x)
	view = min(view, nevents(x))
    
	cat(paste('Dataset: n=', nsamples(x), ', m=', nevents(x), ', |G|=', ngenes(x), '.\n', sep=''))
	cat(paste('Events (types): ', paste(as.types(x), collapse=', '), '.\n', sep=''))
	cat(paste('Colors (plot): ', paste(as.colors(x), collapse=', '), '.\n', sep=''))

	if(has.stages(x))
	{
		cat(paste('Stages: '))
		
		s = unlist(sort(unique(as.stages(x)[, 1])))
		cat((paste(paste(s, collapse=', ', sep=''), '.\n', sep='')))
	}

	cat(paste('Events (', view, ' shown):\n', sep=''))
	cat(paste(paste('\t', rownames(as.events(x)[1: view,]), ':', as.events(x)[1: view, 1], as.events(x)[1: view, 2], sep=' '), collapse='\n'))

	cat(paste('\nGenotypes (', view, ' shown):\n', sep=''))
	head(x$genotypes[,1:view])
}


# Return the number of types in the cohort. 
#
# @x: the dataset.
ntypes = function(x)
{
  return(length(as.types(x)))
}

# Return the number of samples in the cohort. 
#
# @x: the dataset.
nsamples = function(x)
{
	return(nrow(x$genotypes))
}

# Return the number of events in the dataset involving a certain gene or type of event.
#
# @x: the dataset.
# @genes: a list of genes to consider, if NA all genes are used.
# @types: the types of events to consider, if NA all available types are used.
nevents = function(x, genes=NA, types=NA)
{
	return(nrow(as.events(x, genes, types)))
}

# Return the number of genes in the dataset involving a certain type of event.
#
# @x: the dataset.
# @types: the types of events to consider, if NA all available types are used.
ngenes = function(x, types=NA)
{
	return(length(as.genes(x, types=types)))
}

#' Return the number of types in the dataset involving a certain set of genes
#'
#' @param x the dataset.
#' @param genes the genes to consider, if NA all available ones are used.
ntypes = function(x, genes=NA)
{
  return(length(unique(as.events(x, genes=genes)[, 'type'])))
}

# Convert the internal reprensentation of genotypes to numeric, if not. 
#
# @x: the dataset.
enforce.numeric = function(x)
{
  if(!is.numeric(x$genotypes[1,1]))
  {
    rn = as.samples(x)
    x$genotypes = apply(x$genotypes, 2, as.numeric)
    rownames(x$genotypes) = rn
  }
  
  return(x)
}

# Convert the internal reprensentation of genotypes to characters, if not. 
#
# @x: the dataset.
enforce.string = function(x)
{
  if(!is.character(x$genotypes[1,1]))
  {
    rn = as.samples(x)
    x$genotypes = apply(x$genotypes, 2, as.character)
    rownames(x$genotypes) = rn
  }

  return(x)
}

