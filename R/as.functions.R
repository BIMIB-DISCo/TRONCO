# Return all genotypes. Colnames are the keys to access the database
#
#' @title as genotypes
#' @param x: the dataset.
#' @param genes: a list of genes to query, if NA all genes are used.
#' @export
as.genotypes = function(x)
{
	return(x$genotypes)
}

# Return all samples IDs in the cohort.
#
# @x: the dataset.
# @genes: a list of genes to query, if NA all genes are used.
#' @export
as.samples = function(x)
{
	return(rownames(x$genotypes))
}

# Return all gene symbols for which a certain type of event exists.
#
# @x: the dataset.
# @types: the types of events to consider, if NA all available types are used.
#' @export
as.genes = function(x, types=NA)
{
	return(unlist(unique(as.events(x, types=types)[, 'event'])))
}

# Return all events involving certain genes and types.
#
# @x: the dataset.
# @genes: a list of genes to query, if NA all genes are used.
# @types: the types of events to consider, if NA all available types are used.
#' @export
as.events = function(x, genes=NA, types=NA)
{
  ann = x$annotations[, c('type', 'event'), drop=FALSE]

  if(!any(is.na(genes))) ann = ann[ which(ann[, 'event', drop=FALSE] %in% genes) , , drop=FALSE] 
  if(!any(is.na(types))) ann = ann[ which(ann[, 'type', drop=FALSE] %in% types) , , drop=FALSE]   

  return(ann)
}

# Return all stages per sample, if present.
#
# @x: the dataset.
#' @export
as.stages = function(x)
{
  if(has.stages(x)) return(x$stages)
  else return(NA)
}

# Return the (labels) of events in the dataset for a set of genes.
#
# @x: the dataset.
# @genes: a list of genes to query, if NA all available genes are used.
#' @export
as.types = function(x, genes=NA)
{    
  return(unlist(unique(as.events(x, genes=genes)[, 'type'])))
}

# Return all colors map for the events in the dataset.
#
# @x: the dataset.
#' @export
as.colors = function(x)
{
  return(x$types[, 'color'])
}

# Return genotypes for a set of genes and events. 
#
# @x: the dataset.
# @genes: a list of genes to query.
# @types: the types of events to consider, if NA all available types are used.
#' @export
as.gene = function(x, genes, types=NA)
{
  keys = as.events(x, genes=genes, types=types)
  
	data = data.frame(x$genotypes[, rownames(keys)], row.names = as.samples(x))
  colnames(data) = keys[, 'type']
	
	return(data)
}

# For an input dataset merge all the events in a new one.
#
# @x: the dataset.
# @new.type: label for the new type to create
# @new.color: color for the new type to create
#' @export
as.alterations = function(x, new.type = 'Alteration', new.color = 'khaki') {
  merge.types(x, NULL, new.type = new.type, new.color = new.color)
}

# Return true if stages are present
#
# @x: the dataset.
#' @export
has.stages = function(x)
{
  return(! (is.null(x$stages) || is.na(x$stages)))
}


# Return true if duplicate events are present.
#
# @x: the dataset.
#' @export
has.duplicates = function(x) {

  # find duplicate over the dataset
  dup = duplicated(as.events(x))
   
  # return true if at least one duplicate is found
  return(any(dup))
}

# Return duplicated events, if any.
#
# @x: the dataset.
#' @export
duplicates = function(x) {
  as.events(x)[duplicated(as.events(x)),]
}


# Print a short report of a dataset. 
#
# @x: the dataset.
# @view: the number of events to show.
#' @export
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
	to.show = paste( '\t',
				rownames(as.events(x)[1: view,]), ':', 
				as.events(x)[1: view, 1], as.events(x)[1: view, 2], sep=' ')

	cat(paste(to.show, collapse = '\n'))

	cat(paste('\nGenotypes (', view, ' shown):\n', sep=''))
	print(head(x$genotypes[,1:view, drop=FALSE]))
}


# Return the number of types in the cohort. 
#
# @x: the dataset.
#' @export
ntypes = function(x)
{
  return(length(as.types(x)))
}

# Return the number of samples in the cohort. 
#
# @x: the dataset.
#' @export
nsamples = function(x)
{
	return(nrow(x$genotypes))
}

# Return the number of events in the dataset involving a certain gene or type of event.
#
# @x: the dataset.
# @genes: a list of genes to consider, if NA all genes are used.
# @types: the types of events to consider, if NA all available types are used.
#' @export
nevents = function(x, genes=NA, types=NA)
{
	return(nrow(as.events(x, genes, types)))
}

# Return the number of genes in the dataset involving a certain type of event.
#
# @x: the dataset.
# @types: the types of events to consider, if NA all available types are used.
#' @export
ngenes = function(x, types=NA)
{
	return(length(as.genes(x, types=types)))
}

#' Return the number of types in the dataset involving a certain set of genes
#'
#' @param x the dataset.
#' @param genes the genes to consider, if NA all available ones are used.
#' @export
ntypes = function(x, genes=NA)
{
  return(length(unique(as.events(x, genes=genes)[, 'type'])))
}

# Convert the internal reprensentation of genotypes to numeric, if not. 
#
# @x: the dataset.
#' @export
enforce.numeric = function(x)
{
  if(!all(is.numeric(x$genotypes[1,])))
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
#' @export
enforce.string = function(x)
{
  if(!all(is.character(x$genotypes[1,])))
  {
    rn = as.samples(x)
    x$genotypes = apply(x$genotypes, 2, as.character)
    rownames(x$genotypes) = rn
  }

  return(x)
}

# as.pathway - Given a cohort and a pathway, return the cohort with events restricted to genes 
# involved in the pathway. This might contain a new 'pathway' genotype with an alteration mark if
# any of the involved genes are altered.
#
# @pathway.genes: gene - symbols - involved in the pathway
# @pathway.name: pathway name
# @pathway.color: pathway color (for visualization)
# @aggregate.pathway: if TRUE does not show the genes in the pathway
#' @export
as.pathway <- function(x, pathway.genes, pathway.name, 
                       pathway.color='yellow', aggregate.pathway = TRUE) 
{
  is.compliant(x, 'as.pathway: input')
  
  data = x$genotypes
  
  cat(paste('*** Extracting events for pathway: ', pathway.name,'.\n', sep=''))
  
  # Select only those events involving a gene in pathway.genes which is also in x
  y = events.selection(x, NA, filter.in.names=pathway.genes, NA)
  
  # Extend genotypes
  y = enforce.numeric(y)
  
  pathway = data.frame(rowSums(as.genotypes(y)), row.names = as.samples(y), stringsAsFactors = FALSE)  
  pathway[pathway > 1, ] =  1
  colnames(pathway) = pathway.name 

  pathway = import.genotypes(pathway, event.type = 'Pathway', color = pathway.color)

  cat('Pathway extracted succesfully.\n')
  
  if(!aggregate.pathway) pathway = ebind(pathway, y)
  
  if(has.stages(y)) pathway = annotate.stages(pathway, as.stages(y))
  
  is.compliant(pathway, 'as.pathway: output')
  
  return(pathway)
}

sort.by.frequency = function(x)
{
  is.compliant(x)
  
  x = enforce.numeric(x)
  sums = colSums(x$genotypes)
    
  x$genotypes = x$genotypes[, order(sums, decreasing = TRUE), DROP = F]
  x$annotations = x$annotations[colnames(x$genotypes), , DROP = F]

  return(x)  
}

#' Return the number of hypotheses in the dataset
#'
#' @param x the dataset.
#' @export
nhypotheses = function(x)
{
  is.compliant(x$data)
  if ('hstructure' %in% names(x$data$hypotheses)) {
    return(length(ls(x$data$hypotheses$hstructure)))
  }
  return(0)
}

#' Return the name of hypotheses in the dataset
#'
#' @param x the dataset.
#' @export
as.hypotheses = function(x)
{
  is.compliant(x$data)
  if ('hstructure' %in% names(x$data$hypotheses)) {
    return(ls(x$data$hypotheses$hstructure))
  }
}

#' Return all events involving certain genes and types.
#'
#' @param x: the dataset.
#' @param hypotheses: 
#' @export
as.events.hypotheses = function(x, hypotheses=NULL)
{
  is.compliant(x$data)
  ann = x$data$annotations[, c('type', 'event'), drop=FALSE]
  if (is.null(hypotheses)) {
    hypotheses = as.hypotheses(x)
  }
  
  genes_list = NULL
  for(h in hypotheses) {
    if(!h %in% as.hypotheses(x)) {
      stop('Hypothesis ', h, ' not in as.hypotheses(x)')
    }
    g = lapply(colnames(x$data$hypotheses$hstructure[[h]]), function(x){  if(length(i <- grep('^G([0-9]+)$', x))){x[i]}})
    genes_list = append(genes_list, g)
  }  
  genes_list = unique(unlist(genes_list))

  if(!(is.null(genes_list))) ann = ann[ which(rownames(ann) %in% genes_list) , , drop=FALSE] 

  return(ann)
}

# Return all gene symbols for given an hypotheses list
#
# @x: the dataset.
# @hypotheses: the hypotheses to consider, if NULL all available hypotheses are used.
#' @export
as.genes.hypotheses = function(x, hypotheses=NULL) {
  events = as.events.hypotheses(x, hypotheses)
  genes = unique(events[,'event'])
  return(genes)
}
