######## Export to NBS
# Create a .mat file which can be used with NBS clustering
# (ref: http://chianti.ucsd.edu/~mhofree/wordpress/?page_id=26)
#' @param x The TRONCO standard object
#' @param map_hugo_entrez   Hugo_Symbol  - Entrez_Gene_Id
#' @param filename
#' @param filepath 
export.nbs.input = function(x, 
                      map_hugo_entrez,
                      filename = 'tronco_to_nbs.mat',
                      filepath = './')
{
  if (!require(R.matlab)) {
    install.packages('R.matlab', dependencies = TRUE)
    library(R.matlab)
  }
  
  is.compliant(x);
  
  # gene_indiv_mat <- the matrix
  gene_indiv_mat = as.matrix(data)
  
  # remove colnames and rownames from gene_indiv_mat
  rownames(gene_indiv_mat) = NULL
  colnames(gene_indiv_mat) = NULL
  
  # sample_id <- patient id
  sample_id = as.samples(x)
    
  # gene_id_symbol <- sorted name of events
  gene_id_symbol = as.genes(x)
  
  # gene_id_all <- event id positions
  gene_id_all = mapply(function(x) as.numeric(map_hugo_entrez[[which(map_hugo_entrez[,'Hugo_Symbol'] == x), 'Entrez_Gene_Id']]), gene_id_symbol)
  
  filename = if(grepl("\\.mat$", filename)) filename else paste0(filename, ".mat")
  filepath = if(grepl("\\/$", filepath)) filepath else paste0(filepath, "/")
  con = paste0(filepath, filename)
  
  writeMat(con, gene_indiv_mat = gene_indiv_mat, gene_id_all = gene_id_all, sample_id = sample_id, gene_id_symbol = gene_id_symbol)
}

######## Import NBS groups -- current NBS version is XXX

# TODO
