
"cluster.sensitivity" <- function(cluster.map, reference, stages=NA, file=NA) {
  require(gplots)
  require(colorspace)
  require(pheatmap)
  require(RColorBrewer)
  
  # We should use something like this - which does not work
  #requirements = c('gplots', 'colorspace', 'pheatmap', 'RColorBrewer')
  #   aux.fun = function(x) {
  # 	if (!require(x)) {
  # 		cat(paste('Installing required ', x, ' library, with dependencies.\n'))
  #     	install.packages(x, dependencies = TRUE)
  #     	library(x)
  #   	}
  #   	}
  #    lapply(requirents, aux.fun)
  
  if(ncol(cluster.map) == 1) stop('No clustering stability for a unique clustering map!')

  if(!reference %in% colnames(cluster.map)) 
    stop(paste0('The reference cluster specified is not any of: ', 
                paste(colnames(cluster.map), collapse=', '), '.'))

  ref.clust = which(reference == colnames(cluster.map))  
  colnames(cluster.map)[ref.clust] = paste(colnames(cluster.map)[ref.clust],' [reference]',sep='') 
  
  # Transpose data
  cluster.map = t(cluster.map)
  
  # Sort data according to reference row
 cluster.map = cluster.map[, sort(cluster.map[ref.clust, ], decreasing=FALSE, index.return=TRUE)$ix];
  
  # Get unique clustering IDs
  id = apply(cluster.map, 1, unique)
  id = unique(unlist(id))

  # Compute the clustering score
  subdata = cluster.map[-ref.clust,]
  refcol = cluster.map[ref.clust,]
  urefcol = unique(refcol)

  score = rep(0, ncol(cluster.map))
 
  for(i in 1:length(urefcol))
  {
    tmp = as.matrix(subdata[,which(refcol==i)]);
    
    curr.score = 0;
    for (j in 1:nrow(tmp)) 
    {
      curr.cardinality = sort(table(tmp[j,]),decreasing=TRUE)[[1]];
      curr.score = curr.score + (ncol(tmp) - curr.cardinality)/ncol(tmp);			
    }
    
    score[which(refcol==i)] = 1 - (curr.score/nrow(tmp))
  }
   
  # Create annotations
  cn = colnames(cluster.map)

  annotation = data.frame(sensitivity=score, row.names=cn, stringsAsFactors=FALSE)

  if(!all(is.na(stages)))
  	annotation$stage = stages[cn,1]
  	
  # Create colors 
  col = brewer.pal(n = length(id), name = 'Set1')

  different.stages = sort(unique(annotation$stage))
  num.stages = length(different.stages)
  stage.color = append(brewer.pal(n = num.stages, name = 'YlOrRd'), '#FFFFFF')
  names(stage.color) = append(levels(as.factor(different.stages)), NA)

  score.color = brewer.pal(n = 10, name = 'Blues')

  # Annotation colors
  annotation_colors = list(stage=stage.color, sensitivity=score.color)
  
  # Settings
  main = paste0("Sensitivity of clustering assigment for ", reference,', with respect to clusters detected for ', paste(rownames(subdata), collapse=', ')) 
  
  fontsize_col = 3

  pheatmap(cluster.map, 
           scale = "none", 
           cluster_col= F, 
           cluster_rows = F, 
           col=col, 
           main=main,
           fontsize=7,
           fontsize_col=fontsize_col,
           annotation = annotation,
           annotation_colors = annotation_colors,	
           border_color='lightgray',
           border=T,
           margins=c(10,10),
           cellwidth = 3, 
           cellheight = 25,
           #legend_labels = unlist(unique(id)),
           legend=F,
           #	legend_breaks=-1:4
           #treeheight_column = 20,
           filename=file
  )
  
  return(cluster.map)
}
