#### data.edit.R
####
#### TRONCO: a tool for TRanslational ONCOlogy
####
#### See the files COPYING and LICENSE for copyright and licensing
#### information.

rename.type <- function(x, old.name, new.name) {
  # if is compliant x
  is.compliant(x)
  
  if (old.name %in% as.types(x)) {
    x$annotations[ which(x$annotations[,'type'] == old.name), 'type' ] = new.name
    rownames(x$types)[which(rownames(x$types) == old.name)] = new.name
  } else {
    stop(paste(old.name, 'not in as.types(x)'))
  }
  
  is.compliant(x)
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
  print(as.types(x))
  if (type %in% as.types(x)) {
    print('babana')
    drops = rownames(x$annotations[ which(x$annotations[,'type'] == type), ])
    print(drops)
    x$genotypes = subset(x$genotypes, select= -c(unlist(drops)))
    #x$annotations = x$annotations[! (rownames(x$annotations) %in% drops), ]
  } else {
    stop(paste(old.name, 'not in as.types(x)'))
  }
  
  is.compliant(x)
  #return(x)
  
}