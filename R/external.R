##################################################################################
#                                                                                #
# TRONCO: a tool for TRanslational ONCOlogy                                      #
#                                                                                #
##################################################################################
# Copyright (c) 2015, Marco Antoniotti, Giulio Caravagna, Luca De Sano,          #
# Alex Graudenzi, Ilya Korsunsky, Mattia Longoni, Loes Olde Loohuis,             #
# Giancarlo Mauri, Bud Mishra and Daniele Ramazzotti.                            #
#                                                                                #
# All rights reserved. This program and the accompanying materials               #
# are made available under the terms of the GNU GPL v3.0                         #
# which accompanies this distribution                                            #
#                                                                                #
##################################################################################

# Create a .mat file which can be used with NBS clustering
# (ref: https://code.google.com/p/mutex/ )
#' @title export mutex
#' @param x The TRONCO standard object
#' @param filename TODO
#' @param filepath  TODO
#' @param label.mutation TODO
#' @param label.amplification TODO
#' @param label.deletion TODO
#' @export
export.mutex = function(x, 
                      filename = 'tronco_to_mutex',
                      filepath = './',
                      label.mutation = 'SNV',
                      label.amplification = list('High-level Gain'),
                      label.deletion = list('Homozygous Loss'))
{
  
  is.compliant(x)
  data = x
  alteration = list(unlist(label.mutation), unlist(label.amplification), unlist(label.deletion))
  #print(alteration)
  
  # merge amplification
  if (length(label.amplification) >= 0) {
    amplification = label.amplification[[1]]
  }
  if (length(label.amplification) >= 2) {
    amplification = 'amplification'
    data = union.types(data, label.amplification[[1]], label.amplification[[2]], 'amplification', 'red')
  }
  if (length(label.amplification) > 2) {
    for (label in label.amplification[3:length(label.amplification)]) {
      data = union.types(data, label, 'amplification', 'amplification', 'red')
    }
  }
    
  # merge deletion
  if (length(label.deletion) >= 0) {
    deletion = label.deletion[[1]]
  }
  if (length(label.deletion) >= 2) {
    deletion = 'deletion'
    data = union.types(data, label.deletion[[1]], label.deletion[[2]], 'deletion', 'blue')
  }
  if (length(label.deletion) > 2) {
    for (label in label.deletion[3:length(label.deletion)]) {
      data = union.types(data, label, 'deletion', 'deletion', 'blue')
    }
  }
  
  # merge mutation
  if (length(label.mutation) >= 0) {
    mutation = label.mutation[[1]]
  }
  if (length(label.mutation) >= 2) {
    mutation = 'mutation'
    data = union.types(data, label.mutation[[1]], label.mutation[[2]], 'mutation', 'green')
  }
  if (length(label.mutation) > 2) {
    for (label in label.mutation[3:length(label.mutation)]) {
      data = union.types(data, label, 'mutation', 'mutation', 'green')
    }
  }

  samples = rownames(data$genotypes) 
  genes = unique(data$annotation[,'event'])
  
  mutex.matrix = matrix(0, nrow = length(genes), ncol = length(samples))
  colnames(mutex.matrix) = samples
  rownames(mutex.matrix) = genes
  
  #print(mutex.matrix)
  
  # legend:
  # 0: no alteration
  # 1: mutation
  # 2: amplification
  # 4: deletion
  # 3: 1+2 a+m
  # 5: 1+4 d+m
  
  legend = list(1, 2, 4)
  names(legend) = list(mutation, amplification, deletion)
  tronco.matrix = data$genotypes
  
  for(sample in rownames(tronco.matrix)) {
    for(gene in colnames(tronco.matrix)) {
      type = data$annotations[[gene, 'type']]
      
      # print(typeof())
      
      if(type %in% alteration && tronco.matrix[sample, gene] == 1) {
        # print(paste('sample: ', sample, 'event: ', gene, 'type: ', data$annotations[gene, 'type'], 'gene: ', data$annotations[gene, 'event']))
        to.add = legend[[data$annotations[[gene, 'type']]]]
        #print(to.add)
        actual.value = mutex.matrix[data$annotations[[gene, 'event']], sample]
        #print(actual.value)
        mutex.matrix[data$annotations[[gene, 'event']], sample] = actual.value + to.add
      }
    }
  }
  
  # reassign value according to mutex notation
  
  # legend:
  # 0: no alteration
  # 1: mutation
  # 2: amplification
  # 3: deletion
  # 4: 1+2 a+m
  # 5: 1+4 d+m
  
  # move a+m to 10
  mutex.matrix[which(mutex.matrix == 3)] = 10
  
  # move deletion to 3
  mutex.matrix[which(mutex.matrix == 4)] = 3

  # move a+m to 4
  mutex.matrix[which(mutex.matrix == 10)] = 4
    
  mutex.header = append("Symbol", samples)
    
	filepath = if(grepl("\\/$", filepath)) filepath else paste0(filepath, "/")
  con = paste0(filepath, filename)
  write(mutex.header, file=con, sep = "\t", ncolumns = length(mutex.header))
	write.table(mutex.matrix, con, sep="\t", append=T, col.names = F, quote = F)
  
  return(mutex.matrix)
 }




######## Export to NBS
# Create a .mat file which can be used with NBS clustering
# (ref: http://chianti.ucsd.edu/~mhofree/wordpress/?page_id=26)
#' @title export nbs input
#' @param x The TRONCO standard object
#' @param map_hugo_entrez   Hugo_Symbol  - Entrez_Gene_Id
#' @param file TODO
#' @export
export.nbs.input = function(x, 
                      map_hugo_entrez,
                      file = 'tronco_to_nbs.mat')
{

  library(R.matlab)
  
  is.compliant(x);
  
  cat('*** Exporting for NBS v. 0.2\n')

  cat('Preparing binary input matrix\n')
  # gene_indiv_mat <- the matrix
  gene_indiv_mat = as.matrix(x$genotypes)
  
  # remove colnames and rownames from gene_indiv_mat
  rownames(gene_indiv_mat) = NULL
  colnames(gene_indiv_mat) = NULL
  
  cat('Preparing samples IDs \n')
  # sample_id <- patient id
  sample_id = as.samples(x)
    
  cat('Preparing genes list (should be Hugo_Symbol) \n')
  # gene_id_symbol <- sorted name of events
  gene_id_symbol = as.genes(x)
  
  cat('Preparing genes map (should be Hugo_Symbol -> Entrez_Gene_Id) \n')
  if(!('Hugo_Symbol' %in% colnames(map_hugo_entrez))) stop('No Hugo_Symbol column in the input map: ', colnames(map_hugo_entrez))
  if(!('Entrez_Gene_Id' %in% colnames(map_hugo_entrez))) stop('No Entrez_Gene_Id column in the input map: ', colnames(map_hugo_entrez))
  
  gene_id_all = mapply(function(x) as.numeric(map_hugo_entrez[[which(map_hugo_entrez[,'Hugo_Symbol'] == x), 'Entrez_Gene_Id']]), gene_id_symbol)
  
  file = if(grepl("\\.mat$", file)) file else paste0(file, ".mat")
  con = paste0(file)
  
  cat('Writing Matlab file to disk:', file,  ' ..... ' )
  writeMat(con, gene_indiv_mat = gene_indiv_mat, gene_id_all = gene_id_all, sample_id = sample_id, gene_id_symbol = gene_id_symbol)
  cat('DONE')
}

######## Import NBS groups -- current NBS version is XXX

# TODO


######## Export to mutex

## TODO


######## Import mutex groups -- current Mutex version is XXX
# Create a list of unique Mutex groups for a fdr cutoff
# (ref: http://.....)
#' @title import mutex group
#' @param file Mutex results ("ranked-groups.txt" file)
#' @param fdr cutoff for fdr
#' @param display print summary table of extracted groups
#' @export
import.mutex.groups = function(file, fdr=.2, display = TRUE)
{
  # Found somewhere on the web - makes sense
  read.irregular <- function(filenm) 
  {
    fileID <- file(filenm,open="rt")
    nFields <- count.fields(fileID)
    mat <- matrix(nrow=length(nFields),ncol=max(nFields))
    invisible(seek(fileID,where=0,origin="start",rw="read"))
    for(i in 1:nrow(mat) ) {
      mat[i,1:nFields[i]] <-scan(fileID,what="",nlines=1,quiet=TRUE)
    }
    close(fileID)
    df <- data.frame(mat, stringsAsFactors=FALSE)
    return(df)
  }
  
  x = read.irregular(file)
  
  # Check header
  if(any(x[1,1:3] != c('Score', 'q-val', 'Members')))
    warning('File header does not seem to contain \'Score\', \'q-val\' and \'Members field\' - are you
            sure this is a Mutex result file?' )

  # Remove header
  cat(paste('*** Groups extracted - ', (nrow(x) -1), ' total groups.\n', sep=''))
  x = x[-1, , drop = F] # this is c('Score', 'q-val', 'Members')
  x[, 1] = as.numeric(x[,1]) # fdr
  x[, 2] = as.numeric(x[,2]) # q-value
  
  # remove groups  with low fdr
  res = x[which(x[,1] < fdr), , drop = F] 
  
  # remove duplicated groups (permutations)
  res.g = res[, 3:ncol(res)]
  
  for(i in 1:nrow(res.g)) res[i,3:ncol(res)] = sort(res.g[i,], na.last = T)   
    res = res[!duplicated((res[,3:ncol(res), drop=FALSE])), ] 
  
  cat(paste('Selected ', nrow(res), ' unique groups with fdr < ', fdr, '\n', sep=''))
  
  # Create groups
  groups = function(g) {
    # print(g)
    
    # print(g[3:length(g)])
    
    g = g[3:length(g)]
    g = g[!is.na(g)]
    names(g) = NULL
    
    # print(g)
    return(sort(g))
  }
  
  # print(res)
  # print(apply(res, 1, groups))
  
  # G = apply(res, 1, groups)
  # print(G)
  
  # print(res)
  G = list()
  for(i in 1:nrow(res))
  {
    gr = list(groups(res[i, ]))
    names(gr) = paste('MUTEX_GROUP', i, sep='')
    G = append(G,gr)  
  }
  # print(G)
  
  # if(!is.list(G)) 
  # {
    # G = list(as.vector(G))  
    # names(G) = paste('MUTEX_GROUP', 1, sep='')
  # }

  # names(G) = paste('MUTEX_GROUP', 1:length(names(G)), sep='')
  rownames(res) = names(G)
  colnames(res)[1:2] = c('fdr', 'score')
  
  # Summary report
  if(display) 
  { 
    print(res)
  # print(G)
  }
  return(G)
}

#' @export
TCGA.multiple.samples = function(x)
{
  is.compliant(x)
  
  samples = as.samples(x)
  samples.truncated = substring(samples, 0, 12)
  
  patients = unique(samples.truncated)
  
  if(length(patients) != nsamples(x))
  {
    dup.samples.start = which(duplicated(samples.truncated)) 
    dup.samples.last = which(duplicated(samples.truncated, fromLast = T))
    
    return(sort(samples[c(dup.samples.start, dup.samples.last)])) 
  }
  else return(NA)
}

#' @export
TCGA.remove.multiple.samples = function(x)
{
  is.compliant(x, err.fun='Removing TCGA multiple samples (input)')

  dup = TCGA.multiple.samples(x)
  dup.truncated = substring(dup, 0, 12)
  patients = unique(dup.truncated)

  for(i in 1:length(patients))
  {
    patients.samples = which(dup.truncated == patients[i])
    multiple.samples = dup[patients.samples]
    
    cat('Patient', patients[i], 'with sample aliquotes\n' )
    print(substring(multiple.samples, 14, 29))

    keep = max(multiple.samples)
    discard = multiple.samples[which(multiple.samples != keep)]
    
    cat('Selecting', keep, '\n')
    x = delete.samples(x, discard)
  }
    
  is.compliant(x, err.fun='Removing TCGA multiple samples (output)')
  return(x)
}

#' @export
TCGA.shorten.barcodes = function(x)
{
  is.compliant(x, err.fun='Shartening TCGA barcodes (input)')

  # Check if it has duplicated barcodes
  if(!all(is.na(TCGA.multiple.samples(x))))
    stop(
      paste('This dataset contains multiple samples for some patients - cannot consolidate.',
            '\n Samples with barcodes indicating multiple patients: \n', paste(TCGA.multiple.samples(x), collapse = '\n'), '.'
            , sep =''))
  
  # Shorten sample barcodes
  rownames(x$genotypes) = substring(rownames(x$genotypes), 0, 12)
  if(has.stages(x)) rownames(x$stages) = rownames(x$genotypes)

  is.compliant(x, err.fun='Shartening TCGA barcodes (output)')
  return(x)    
}
  
#' @export
TCGA.map.clinical.data = function(file, sep='\t', column.samples, column.map)
{
  
  data = read.delim(
    file = file,
    sep = sep,
    header = TRUE,
    stringsAsFactors=F)
  
  if(!(column.samples %in% colnames(data))) 
    stop(paste('Cannot find samples column \"', column.samples, '\". Available columns: \n\t',
               paste(colnames(data), collapse='\n\t'), sep = ''))
  
  if(!(column.map %in% colnames(data))) 
    stop(paste('Cannot find required map column \"', column.map, '\". Available columns: \n\t',
               paste(colnames(data), collapse='\n\t'), sep = ''))
  
  map = data.frame(data[, column.map], row.names = data[, column.samples])
  colnames(map) = column.map
  
  return(map)
}

sample.RColorBrewer.colors = function(palette, ncolors)
{
  if(!palette %in% rownames(brewer.pal.info)) stop('Invalid RColorBrewer palette.')

  pmax.cols = brewer.pal.info[palette, 'maxcolors']
  
  cols = min(pmax.cols , ncolors)
  cols = ifelse(cols < 3, 3, cols)
  
  colors = brewer.pal(n=cols, name=palette)
  if(ncolors < 3) colors = colors[1:ncolors]
  else colors =  colorRampPalette(colors)(ncolors)
  
  return(colors)
}
