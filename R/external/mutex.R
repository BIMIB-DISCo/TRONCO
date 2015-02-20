######## Export to mutex
# Create a .mat file which can be used with NBS clustering
# (ref: http://chianti.ucsd.edu/~mhofree/wordpress/?page_id=26)
#' @param x The TRONCO standard object
#' @param map_hugo_entrez   Hugo_Symbol  - Entrez_Gene_Id
#' @param filename
#' @param filepath 
export.nbs = function(x, 
                      map_hugo_entrez,
                      filename = 'tronco_to_nbs.mat',
                      filepath = './')
{
  if (!require(R.matlab)) {
    install.packages('R.matlab', dependencies = TRUE)
    library(R.matlab)
  }
  
  is.compliant(x);
  
  data = x$genotypes
  event = x$annotation[, 'event']
  # gene_indiv_all <- the matrix
  gene_indiv_all = as.matrix(data)
  
  # sample_id <- patient id
  sample_id = rownames(gene_indiv_all)
  
  # remove colnames and rownames from gene_indiv_mat
  rownames(gene_indiv_all) = NULL
  colnames(gene_indiv_all) = NULL
  
  # gene_id_symbol <- sorted name of events
  gene_id_symbol = mapply(function(x) event[[x]], colnames(data))
  
  # gene_id_all <- event id positions
  gene_id_all = mapply(function(x) map_hugo_entrez[[which(map_hugo_entrez[,'Hugo_Symbol'] == x), 'Entrez_Gene_Id']], gene_id_symbol)
  
  filename = if(grepl("\\.mat$", filename)) filename else paste0(filename, ".mat")
  filepath = if(grepl("\\/$", filepath)) filepath else paste0(filepath, "/")
  con = paste0(filepath, filename)
  
  writeMat(con, gene_indiv_all = gene_indiv_all, gene_id_all = gene_id_all, sample_id = sample_id, gene_id_symbol = gene_id_symbol)
}

######## Import mutex groups -- current Mutex version is XXX
# Create a list of unique Mutex groups for a fdr cutoff
# (ref: http://.....)
#' @param file Mutex results ("ranked-groups.txt" file)
#' @param fdr cutoff for fdr
#' @param display print summary table of extracted groups
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
  x = x[-1, ] # this is c('Score', 'q-val', 'Members')
  x[, 1] = as.numeric(x[,1]) # fdr
  x[, 2] = as.numeric(x[,2]) # q-value
  
  # remove groups  with low fdr
  res = x[which(x[,1] < fdr), ] 
  
  # remove duplicated groups (permutations)
  res.g = res[, 3:ncol(res)]
  for(i in 1:nrow(res.g)) res[i,3:ncol(res)] = sort(res.g[i,], na.last = T)   
    res = res[!duplicated((res[,3:ncol(res)])), ] 
  
  cat(paste('Selected ', nrow(res), ' unique groups with fdr < ', fdr, '\n', sep=''))
  
  # Create groups
  groups = function(g) {
    g = g[3:length(g)]
    g = g[!is.na(g)]
    names(g) = NULL
    
    return(sort(g))
  }
  
  G = apply(res, 1, groups)
  names(G) = paste('MUTEX_GROUP', 1:length(names(G)), sep='')
  rownames(res) = names(G)
  colnames(res)[1:2] = c('fdr', 'q-score')
  
  # Summary report
  if(display) print(res)
 
  return(G)
}
