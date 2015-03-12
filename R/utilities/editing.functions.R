#### data.edit.R
####
#### TRONCO: a tool for TRanslational ONCOlogy
####
#### See the files COPYING and LICENSE for copyright and licensing
#### information.

rename.type <- function(x, old.name, new.name) {
  # if is compliant x
  is.compliant(x, 'rename.type: pre')
  types = as.types(x)
  
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
  print(x$types)
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
    x$types = t(t(x$types[ which(rownames(x$types) %in% unique(x$annotations[,"type"])), ]))
    colnames(x$types) = 'color'
  } else {
    stop(paste(gene, 'not in as.types(x)'))
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
    x$stages = x$stages[!rownames(x$stages) %in% del, ]
  }
  
  is.compliant(x)
  
  x = trim(x)
  return(x)
}
