#### data.edit.R
####
#### TRONCO: a tool for TRanslational ONCOlogy
####
#### See the files COPYING and LICENSE for copyright and licensing
#### information.

consolidate.data = function(x, print = FALSE){
  is.compliant(x)
  ind = list()
  zeros = list()
  ones = list()
  
  for(i in 1:nevents(x)) {
    ev = list()
    
    for(j in i:nevents(x)) {
      if(i != j && all(x$genotypes[, i] == x$genotypes[, j]) && !j %in% ind )
        ev = append(ev, j)
    }
    
    if(length(ev) > 0)
    {      
      ind.pool = rbind(as.events(x)[ c(i, unlist(ev)),])

      if(all(x$genotypes[, i] == 1)) ones = append(ones, list(ind.pool))
      if(all(x$genotypes[, i] == 0)) zeros = append(zeros, list(ind.pool))
      if(sum(x$genotypes[, i] < nsamples(x))) ind = append(ind, list(ind.pool))
        
      if(print){
        
        if(all(x$genotypes[, i] == 1)) cat('\nEvents altered across all samples:\n')
        if(all(x$genotypes[, i] == 0)) cat('\nEvents with no alterations across samples:\n')
        if(sum(x$genotypes[, i] < nsamples(x))) cat('\nIndistinguishable events:\n')
        
        print(ind.pool)
        cat('Total number of events for these genes: ')
        cat(paste(nevents(x, as.events(x, genes=c(ind.pool)))), '\n')
        
#        cat('All events for these genes are:\n')
#        print(as.events(x, genes=c(ind.pool)))	
      }
    }
  }

  

  ret = NULL
  ret$indistinguishable = ind
  ret$zeroes = zeros
  ret$ones = ones

  #names(ind) = paste0('INDISTINGUISHABLE', 1:length(ind))


  return(ret)
}

#' Annotate a descriotion on the selected dataset
#'
#' @param x A TRONCO compliant dataset.
#' @param label A string
#' @return A TRONCO genotypes matrix.
#' @export annotate.description
annotate.description = function(x, label)
{
  if(as.description(x) != "")
    warning(paste('Old description substituted: ', as.description(x), '.'))
  
  x$name = label
  return(x)
}

#' @export
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
  cat('Events of type', old.name, 'renamed as', new.name, '.\n')
  
  is.compliant(x, err.fun = 'rename.type: output')
  return(x)
}

#' @export
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

#' @export
delete.type <- function(x, type) {
  # if is compliant x
  is.compliant(x)

  if(has.model(x)) {
    stop("There's a reconstructed model, a type cannot be deleted now. \nUse delete.model()")
  }

  for(pattern in as.patterns(x)) {
    if(type %in% as.types.in.patterns(x, patterns=pattern)) {
      stop('Found type \"', type, '\" in pattern \"', pattern, '\". Delete that pattern first.\n')
    }
  }
  
  if (type %in% as.types(x)) {
    events = as.events(x, types=setdiff(as.types(x), type))

    x$genotypes = as.genotypes(x)[,rownames(events), drop=F]
    x$annotations = events
    x$types = x$types[which(rownames(x$types) != type), , drop=F]
  } else {
    stop(paste(type, 'not in as.types(x)'))
  }
  
  is.compliant(x)
  return(x)
}

#' @export delete.gene
delete.gene <- function(x, gene) {
  # if is compliant x
  is.compliant(x, 'delete:gene: input')

  if(has.model(x)) {
    stop("There's a reconstructed model, a type cannot be deleted now. \nUse delete.model()")
  }
  
  if (all(gene %in% as.genes(x))) {
    
    for(pattern in as.patterns(x)) {
      for(g in gene) {
        if(g %in% as.genes.in.patterns(x, patterns=pattern)) {
          stop('Found gene \"', g, '\" in pattern \"', pattern, '\". Delete that pattern first.\n')
        }
      }
    }


    drops = rownames(as.events(x, genes = gene))
    x$genotypes = x$genotypes[, -which( colnames(x$genotypes) %in% drops )]
    x$annotations = x$annotations[ -which (rownames(x$annotations) %in% drops), ]
    
    # TO DO: something better than this t(t(...))
    x$types = x$types[ which(rownames(x$types) %in% unique(x$annotations[,"type"])), , drop=FALSE]
    colnames(x$types) = 'color'

  } else {
    stop(paste(gene[!(gene %in% as.genes(x))], collapse= ','),  'are not in the dataset -- as.genes(x).')
  }
  
  is.compliant(x, 'delete:gene: output')
  return(x)
}

#' @export delete.event
delete.event <- function(x, gene, type) {

  for(pattern in as.patterns(x)) {
    events = as.events.in.patterns(x, patterns=pattern)
    if(length(which(events[,'type'] == type & events[,'event'] == gene)) > 0) {
      stop('Found event \"(', gene, ', ', type, ')\" in pattern \"', pattern, '\". Delete that pattern first.\n')
    }
  }
  
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

#' Delete of a single hypothesis
#' @export delete.hypothesis
delete.hypothesis = function(x, event=NULL, cause=NULL, effect=NULL)
{
  if(has.model(x)) {
    stop("There's a reconstructed model, hypotheses cannot be deleted now. \nUse delete.model()")
  }

  hypo_map = as.hypotheses(x)
  to_remove = c()
  if(!is.null(event)) {
    if(length(event) == 1 && event %in% as.events(x)[,'event']) {
    cause_del = which(hypo_map[,'cause event'] == event)
    effect_del = which(hypo_map[,'effect event'] == event)
    to_remove = unique(c(to_remove, cause_del, effect_del))
    } else {
      stop('Select only one event present in as.events')
    }
  }

  if(! is.null(cause)) {
    if(cause %in% as.events(x)[,'event']) {
    cause_del = which(hypo_map[,'cause event'] == cause)
    to_remove = unique(c(to_remove, cause_del))
    } else {
      stop('Wrong cause, select only events present in as.events')
    }
  }

  if(! is.null(effect)){
    if( effect %in% as.events(x)[,'event']) {
    effect_del = which(hypo_map[,'effect event'] == effect)
    to_remove = unique(c(to_remove, effect_del))
    } else {
      stop('Wrong effect, select only events present in as.events')
    }
  }

  x$hypotheses$num.hypotheses = x$hypotheses$num.hypotheses - 1
  x$hypotheses$hlist = x$hypotheses$hlist[-to_remove, ,drop=F]

  is.compliant(x)

  return(x)
}

#' @export delete.pattern
delete.pattern = function(x, pattern) {
  if(has.model(x)) {
    stop("There's a reconstructed model, a pattern cannot be deleted now. \nUse delete.model()")
  }

  if(! pattern %in% as.patterns(x)) {
    stop(paste(pattern, " not in as.patterns()"))
  }

  x = delete.hypothesis(x, pattern)
  x$annotations = x$annotations[-which(rownames(x$annotations) == pattern), , drop=F]
  x$genotypes = x$genotypes[, -which(colnames(x$genotypes) == pattern), drop=F]
  x$hypotheses$patterns[pattern] = NULL
  x$hypotheses$pvalues = NULL

  for(atom in names(x$hypotheses$atoms)) {
    if(x$hypotheses$atoms[atom] == pattern) {
      x$hypotheses$atoms[atom] = NULL
    }
  }

  rm(list=pattern, envir=x$hypotheses$hstructure)

  # chiave da togliere model$hypotheses$`gene x` se contenuto = 'XOR_EZH2'
  # remove o rm?

  is.compliant(x)

  return(x)
}

#' @export delete.model
delete.model = function(x) {
  if (! has.model(x)) {
    stop("No model to delete in dataset")
  }
  x$model = NULL
  x$confidence = NULL
  x$parameters = NULL
  x$adj.matrix.prima.facie = NULL
  x$execution.time = NULL

  is.compliant(x)

  return(x)
}


#' @export
change.color = function(x, type, new.color)
{
  is.compliant(x)
  
  x$types[type, ] = new.color
  
  is.compliant(x)
  return(x)
}

#' @export
delete.samples = function(x, samples) {
  is.compliant(x, 'delete.samples input')
  stages = has.stages(x)
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
  
  x$genotypes = x$genotypes[!rownames(x$genotypes) %in% del, , drop = F ]
  
  if(stages) {
    x$stages = x$stages[!rownames(x$stages) %in% del, , drop=FALSE]
  }
  
  is.compliant(x, 'delete.samples output')
  
  return(x)
}


# intersect.genomes = F -> just samples 
#' @export
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


#' @export annotate.stages
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


# Binds events from one or more datasets, which must be defined over the same set of samples.
# @... the input datasets
#' @export
ebind = function(...)
{
  # merge two  datasets at a time.
  events.pairwise.bind = function(x, y)
  {  
    is.compliant(x, 'ebind: input x')
    is.compliant(y, 'ebind: input y')
  
    samples.intersect = intersect(as.samples(x), as.samples(y))
    if(!(setequal(samples.intersect, as.samples(x)) && setequal(samples.intersect, as.samples(y))))
      stop('Datasets have different samples, won\'t bind!')
    
    z = list()

  y$genotypes = y$genotypes[rownames(x$genotypes), , drop = FALSE]
    y$stages = y$stages[rownames(y$genotypes), , drop=FALSE]
     
    # Copy genotype matrix, and sets its rownames (samples)
    z$genotypes = cbind(x$genotypes, y$genotypes)
    colnames(z$genotypes) = paste('G', 1:ncol(z$genotypes), sep='')
    
    # Copy annotations for gene symbols etc.
    z$annotations = rbind(x$annotations, y$annotations)
    rownames(z$annotations) = colnames(z$genotypes)
    
    # Copy types
    z$types = unique(rbind(x$types, y$types)) 

    # print(x)
    # print(y)
    # print(has.stages(x))
    # print(has.stages(y))
    # print(as.stages(x) == as.stages(y))
    
    # Copy stages, if present 
    if(has.stages(x) && has.stages(y))
    {
        stages.x = as.stages(x)
        stages.y = as.stages(y)
        
        stages.x = stages.x[!is.na(stages.x)]
      stages.y = stages.y[!is.na(stages.y)]
      
      if(any(stages.x != stages.y))
          stop('Patients have different stages, won\'t merge!')
    }
    
    if(has.stages(x)) 
        z = annotate.stages(z, as.stages(x))

    is.compliant(z, 'ebind: output')
    
    return(z)
  }
  
  
  input = list(...)

  cat('*** Binding events for', length(input), 'datasets.\n')
  return(Reduce(events.pairwise.bind, input))
 }

# Binds samples from one or more datasets, which must be defined over the same set of events
# @... the input datasets
#' @export
sbind = function(...)
{
  # merge two datasets at a time.
  samples.pairwise.bind = function(x, y)
  {  
    is.compliant(x, 'sbind: input x')
    is.compliant(y, 'sbind: input y')
    
    if(!all(as.events(x) == as.events(y)))
    {
     # cat(as.events(x)[as.events(x) == as.events(y),])
      stop('Datasets have different events, can not bind!')
    }
    z = list()
    
    # Copy genotypes and annotations
    z$genotypes = rbind(x$genotypes, y$genotypes)
    z$annotations = x$annotations
    z$types = x$types
        
    # Copy stages, if present
    if(has.stages(x) || has.stages(y)) 
    {
      if(has.stages(x)) xstages = as.stages(x)
      else xstages = matrix(rep('NA', nsamples(x)), nrow=nsamples(x))
      
      if(has.stages(y)) ystages = as.stages(y)
      else ystages = matrix(rep('NA', nsamples(y)), nrow=nsamples(y))
      
      z$stages = (rbind(x$stages, y$stages))
      
      colnames(z$stages) = 'stage'
    }
    
    
    is.compliant(z, 'sbind: output')
    
    return(z)
  }
  
  input = list(...)
  return(Reduce(samples.pairwise.bind, input))
}

# For an input dataset merge all the events of two or more distincit types (e.g., say that missense and indel
# mutations are events of a unique "mutation" type)
#
# @ x : the input dataset
# @ ...: type to merge
# @ new.type: label for the new type to create
# @ new.color: color for the new type to create
#' @export merge.types
merge.types = function(x, ..., new.type = "new.type", new.color = "khaki") {

  # #   # internal function, merge two types
  # merge.two.types = function(x, type.one, type.two, new.type = paste(type.one, type.two, sep = ":"), new.color = "khaki") {
# is.compliant(x, "merge.types: input x")

  # cat("Merging events of type", type.one, "(", nevents(x, types = type.one), ") and", type.two, "(", nevents(x, types = type.two), ")-->", new.type, 
  # "\n")

  # # The returned object
  # z = x

  # # Events of type one/two
  # ev.tone = as.events(x, types = type.one)
# ev.ttwo = as.events(x, types = type.two)

  # # We need the list of genes which have events of both type one and two, and which are specific to a certain type
  # shared.genes = intersect(ev.ttwo[, "event"], ev.tone[, "event"])
# tone.genes = setdiff(ev.tone[, "event"], ev.ttwo[, "event"])
# ttwo.genes = setdiff(ev.ttwo[, "event"], ev.tone[, "event"])

  # #   print(shared.genes)
  # #   print(tone.genes)
# #print(ttwo.genes)

  # # Change names of events which are type specific to type one
  # if (length(tone.genes) > 0) {
# keys.tone.genes = rownames(as.events(x, genes = tone.genes, types = type.one))
# z$annotations[keys.tone.genes, "type"] = new.type
# }

  # # Change names of events which are type specific to type two
  # if (length(ttwo.genes) > 0) {
# keys.ttwo.genes = rownames(as.events(x, genes = ttwo.genes, types = type.two))
# z$annotations[keys.ttwo.genes, "type" ] = new.type
# }

  # # Add new.type to the $types annotation if it does not exist - this leaves "z" still compliant
  # if (!(new.type %in% rownames(z$types))) {
# z$types = rbind(z$types, new.color)
# rownames(z$types)[length(z$types)] = new.type
# }

  # # Now, work with genes which have both types of events
  # if (length(shared.genes) > 0) {
# # These are genotypes restricted to the events we want to process
# data.tone = as.gene(x, genes = shared.genes, types = type.one)
# data.ttwo = as.gene(x, genes = shared.genes, types = type.two)

  # # We can delete the genotypes which we just extracted 
  # for (i in 1:length(data.tone)) z = delete.event(z, gene = shared.genes[i], type = type.one)
# for (i in 1:length(data.tone)) z = delete.event(z, gene = shared.genes[i], type = type.two)

  # for (i in 1:ncol(data.tone)) {
  # # Get the 1s for each event         
# one.tone = data.tone[, i, drop = F] == 1
# one.ttwo = data.ttwo[, i, drop = F] == 1

  # # Build a simple matrix for that
  # geno.tone.ttwo = matrix(rep(0, nsamples(x)), ncol = 1)
# geno.tone.ttwo[which(one.tone), ] = 1
# geno.tone.ttwo[which(one.ttwo), ] = 1

  # # With its colnames 
  # colnames(geno.tone.ttwo) = shared.genes[i]
# rownames(geno.tone.ttwo) = as.samples(x)

  # # Create TRONCO input and bind datasets         
  # genotype = import.genotypes(geno.tone.ttwo, event.type = new.type, color = new.color)
# z = ebind(z, genotype)
# }

  # } 
  
  # # remove types we just merged
  # z$types = z$types[unique(z$annotations[, 'type']), , drop = FALSE]

  # is.compliant(z)
  # return(z)
# }

  # check if x is compliant
  is.compliant(x)

  input = list(...)

  # TODO Change this in a better way (deafult ellipsis?)  
  if (length(input) == 1 && is.null(input[[1]])) {
    input = as.list(as.types(x))
  }

  cat(paste("*** Aggregating events of type(s) {", paste(unlist(input), collapse = ", ", sep = ""), "} in a unique event with label \"", new.type, "\".\n", sep = ""))


  if (length(input) <= 1) {
    cat("One input type provided, using renaming functions.\n")

    x = rename.type(x, input[[1]], new.type)
    x = change.color(x, new.type, new.color)
    return(x)
  }


  types.check = lapply(input, function(type) {
    type %in% as.types(x)
  })

  if (!all(unlist(types.check))) {
    t = (input[!unlist(types.check)])

    stop("No events of type '", t, "' in input dataset, will not merge.")
  }

  if (any(duplicated(input))) {
    stop("No duplicated types are allowed, will not merge.")
  }

  if (new.type %in% as.types(x)) {
    stop(paste0(new.type, "is already used in input dataset, will not merge"))
  }

  input = unlist(input)



  genes = as.genes(x, types = input)
  cat("Dropping event types", paste(input, collapse = ", ", sep = ""), "for", length(genes), "genes.\n")
  geno.matrix = matrix(, nrow = nsamples(x), ncol = length(genes))

  pb = txtProgressBar(1, length(genes), style = 3)
  flush.console()


  for (i in 1:length(genes)) {
    setTxtProgressBar(pb, i)

    geno = as.matrix(rowSums(as.gene(x, genes[i], types = input)))
    geno[geno > 1] = 1

    geno.matrix[, i] = geno
  }


  rownames(geno.matrix) = as.samples(x)
  colnames(geno.matrix) = genes
  close(pb)


  z = import.genotypes(geno.matrix, event.type = new.type, color = new.color)
  if (has.stages(x)) 
    z = annotate.stages(z, as.stages(x))


  y = x
  for(i in input)
  {y = delete.type(y, i)

	# print(show(y))
}

  w = ebind(y, z)
  is.compliant(w)

  return(w)

}

# Deletes all events which have frequency 0 in the dataset.
#' @export
trim = function(x) {
  is.compliant(x, 'trim: input')
  
  x = enforce.numeric(x)
  
  del = names(which(colSums(x$genotypes) == 0))

  x$genotypes = x$genotypes[, !colnames(x$genotypes) %in% del, drop = FALSE]
  x$annotations = x$annotations[!rownames(x$annotations) %in% del, , drop = FALSE]

  x$types = matrix(x$types[ unique(x$annotations[, 'type', drop = FALSE]), , drop = FALSE ], ncol=1)
  rownames(x$types) = unique(x$annotations[,'type', drop = FALSE])
  colnames(x$types) = 'color'
  

  is.compliant(x, 'trim: output')   
  return(x)
}

# TODO: check
# Split cohort (samples) into groups, return either all groups or a specific group.
#
# x: cohort
# clusters: groups, a map  for each sample in 'x' to a group.  
#' @export
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
    
    cat(paste('Group \"', cluster.labels[i,1], '\" has ', length(samples.in.cluster), 
              ' samples.\n', sep=''))
    
    y$genotypes = data[samples.in.cluster, ];
    y$annotations = x$annotations
    y$types = x$types
        
    if(!is.null(x$stages)) {
      y$stages = as.matrix(x$stages[samples.in.cluster, ]) 
      rownames(y$stages) = samples.in.cluster
      }
    
    is.compliant(y, 'subtypes.split partitionig')
    
    
    
    partitions = append(partitions, list(y)) 
    names(partitions)[i] = cluster.labels[i,1]
  }
    
  return(partitions)
}

# Split events into groups according to their types.
#
# x: cohort
#' @export
tsplit <- function(x) 
{
# Parametro per estrarre un tipo solo di evento?
}

merge.events = function(x, events, new.event, new.type, event.color)
{
	if(ncol(events) != 2) stop('ERR - badformed events')
	
	x = enforce.numeric(x)
	
	x.ev = NULL 
	y = x
	for(i in 1:nrow(events))
	{
		x.ev = rbind(x.ev, as.events(x, genes = events[i,2], types = events[i,1]))
		y = delete.event(y, gene = events[i,2], type = events[i,1])
	}
	
	genos = x$genotypes[, rownames(x.ev)]
	genos = matrix(rowSums(genos), nrow = nsamples(x))
	genos[genos > 1] = 1
	colnames(genos) = new.event
	rownames(genos) = as.samples(x)

	genos = import.genotypes(genos, event.type = new.type, color = event.color)
	y = ebind(y, genos)

	# print(as.events(y))
	return(y)
}
