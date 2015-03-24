#### data.edit.R
####
#### TRONCO: a tool for TRanslational ONCOlogy
####
#### See the files COPYING and LICENSE for copyright and licensing
#### information.

rename.type <- function(x, old.name, new.name) {
  is.compliant(x, 'rename.type: input dataset')
  types = as.types(x)
  
  if(old.name == new.name) return(x)	
  
  if (old.name %in% types) {
    x$annotations[ which(x$annotations[,'type'] == old.name), 'type' ] = new.name
    if(! new.name %in% types) {
      rownames(x$types)[which(rownames(x$types) == old.name)] = new.name
    } else {
      x$types = x$types[!rownames(x$types) %in% list(old.name), , drop=F]
    }
  } else {
    stop(paste(old.name, 'not in as.types(x)'))
  }
  cat('Events of type', old.name, 'renamed as', new.name, '')
  
  is.compliant(x, err.fun = 'rename.type: output')
  return(x)
}

rename.gene <- function(x, old.name, new.name) {
  # if is compliant x
  is.compliant(x)
  
  if (old.name %in% as.genes(x)) {
    x$annotations[ which(x$annotations[,'event'] == old.name), 'event' ] = new.name
    
  } else {
    stop(paste(old.name, 'not in as.genes(x)'))
  }
  
  is.compliant(x)
  return(x)
}

delete.type <- function(x, type) {
  # if is compliant x
  is.compliant(x)
  
  if (type %in% as.types(x)) {
    
    drops = rownames(x$annotations[ x$annotations[, "type"] == type,])
    x$genotypes = x$genotypes[, -which( colnames(x$genotypes) %in% drops )]
    x$annotations = x$annotations[ -which (rownames(x$annotations) %in% drops), ]
    
    # TO DO: something better than this t(t(...))
    x$types = matrix(x$types[ which(rownames(x$types) != type), ])
    rownames(x$types) = unique(x$annotations[,'type'])
    colnames(x$types) = c('color')
  } else {
    stop(paste(type, 'not in as.types(x)'))
  }
  
  is.compliant(x)
  return(x)
}

delete.gene <- function(x, gene) {
  # if is compliant x
  is.compliant(x, 'delete:gene: input')
  
  if (gene %in% as.genes(x)) {
    
    drops = rownames(as.events(x, genes = gene))
    x$genotypes = x$genotypes[, -which( colnames(x$genotypes) %in% drops )]
    x$annotations = x$annotations[ -which (rownames(x$annotations) %in% drops), ]
    
    # TO DO: something better than this t(t(...))
    x$types = x$types[ which(rownames(x$types) %in% unique(x$annotations[,"type"])), , drop=FALSE]
    colnames(x$types) = 'color'
  } else {
    stop(paste(gene, 'not in as.types(x)'))
  }
  
  is.compliant(x, 'delete:gene: output')
  return(x)
}

delete.event <- function(x, gene, type) {
  
  is.compliant(x, 'delete.event: input')
  
  if (all(c(type, gene) %in% as.events(x))) 
  {    
    drops = rownames(as.events(x, genes = gene, types = type))
    x$genotypes = x$genotypes[, -which( colnames(x$genotypes) %in% drops ), drop = FALSE]
    x$annotations = x$annotations[ -which (rownames(x$annotations) %in% drops), , drop = FALSE]
    
    # TO DO: something better than this t(t(...))
    x$types = x$types[ which(rownames(x$types) %in% unique(x$annotations[,"type"])), , drop=FALSE]
  } 
  else {
    stop(paste(type, gene, ' - not in as.events(x)'))
  }
  
  is.compliant(x, 'delete:gene: output')
  return(x)
}


change.color = function(x, type, new.color)
{
  is.compliant(x)
  
  x$types[type, ] = new.color
  
  is.compliant(x)
  return(x)
}

delete.samples = function(x, samples) {
  is.compliant(x)
  
  del = list()
  actual.samples = as.samples(x)
  samples = unique(samples)
  for (sample in samples) {    
    if(!sample %in% actual.samples) {
      warning('Sample: ', sample, ' not in as.samples(x)')
    } else {
      del = append(del, sample)
    }
  }
  
  x$genotypes = x$genotypes[!rownames(x$genotypes) %in% del, ]
  
  if("stages" %in% names(x)) {
    x$stages = x$stages[!rownames(x$stages) %in% del, , drop=FALSE]
  }
  
  is.compliant(x)
  
  return(x)
}


# intersect.genomes = F -> just samples 
intersect.datasets = function(x,y, intersect.genomes = TRUE)
{
  is.compliant(x)
  is.compliant(y)
  
  # Common samples and genes (according to intersect.genomes)
  samples = intersect(as.samples(x), as.samples(y))
  genes = ifelse(intersect.genomes, 
                 intersect(as.genes(x), as.genes(y)), #  intersect.genomes -> INTERSECTION
                 unique(c(as.genes(x), as.genes(y)))) # !intersect.genomes -> UNION
  
  report = data.frame(row.names = c('Samples', 'Genes'))
  report$x = c(nsamples(x), ngenes(x))
  report$y = c(nsamples(y), ngenes(y))
  
  # Restrict genes - only if intersect.genomes = T
  if(intersect.genomes)
  {
    x = events.selection(x, filter.in.names=genes)
    y = events.selection(y, filter.in.names=genes)
  }
  
  # TODO: check they have no events in common!
  # if(as.events(x) )
  
  # Restric stamples
  x = samples.selection(x, samples)
  y = samples.selection(y, samples)
  
  # Result
  z = ebind(x,y)
  
  
  cat('*** Intersect dataset [ intersect.genomes =', intersect.genomes, ']\n')
  report$result = c(nsamples(z), ngenes(z))
  print(report)
  
  return(z)    
}



annotate.stages = function(x, stages, match.TCGA.patients = FALSE)
{
  if(is.null(rownames(stages))) stop('Stages have no rownames - will not add annotation.')
  
  samples = as.samples(x)
  
  # Just for temporary - will be shortened to make a simple check...
  samples.temporary = samples
  if(match.TCGA.patients) samples.temporary = substring(samples.temporary, 0, 12)
  
  if(!any(samples.temporary %in% rownames(stages))) 
    stop('There are no stages for samples in input dataset - will not add annotation.')  
  
  # Notify if something gets lost
  if(has.stages(x)) warning('Stages in input dataset overwritten.')
  
  # Actual stages
  x$stages = data.frame(row.names = samples, stringsAsFactors=FALSE)
  x$stages[, 'stage'] = as.character(NA)
  
  for(i in 1:nsamples(x))
  {
    if(!match.TCGA.patients)
      x$stages[i, ] = as.character(stages[as.samples(x)[i], ])
    else 
    {    
      # Potential match if x's samples are long TCGA barcodes and stages are TCGA patients barcdodes (short)
      short.name = substring(as.samples(x)[i], 0 , 12)
      x$stages[i, 'stage'] = as.character(stages[short.name, ])
    }
  }
  
  count.na = is.na(x$stages)
  if(any(count.na)) warning(paste(length(which(count.na)), ' missing stages were added as NA.'))
  
  return(x)
}
