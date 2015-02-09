# Create a .mat file which can be used with NBS clustering
# (ref: https://code.google.com/p/mutex/ )
#' @param x The TRONCO standard object
#' @param filename
#' @param filepath 
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