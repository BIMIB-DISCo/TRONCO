######## Export to NBS
# Create a .mat file which can be used with NBS clustering
# (ref: http://chianti.ucsd.edu/~mhofree/wordpress/?page_id=26)
#' @param x The TRONCO standard object
#' @param map_hugo_entrez   Hugo_Symbol  - Entrez_Gene_Id
#' @param filename
#' @param filepath 
export.nbs.input = function(x, 
                      map_hugo_entrez,
                      file = 'tronco_to_nbs.mat')
{
  if (!require(R.matlab)) {
    install.packages('R.matlab', dependencies = TRUE)
    library(R.matlab)
  }
  
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
