
"cluster.sensitivity" <- function(cluster.map, reference, dataset) {
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

  if(!reference %in% colnames(data)) 
    stop(paste0('The reference cluster specified is not any of: ', 
                paste(colnames(data), collapse=', '), '.'))

  ref.clust = which(reference == colnames(data))  
  colnames(cluster.map)[ref.clust] = paste(colnames(cluster.map)[ref.clust],' [ref.]',sep='') 
  
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

  score = rep(0,ncol(cluster.map))
 
 


  for(i in 1:length(urefcol))
  {
    tmp = as.matrix(subdata[,which(refcol==i)]);
    
    curr.score = 0;
    for (j in 1:nrow(tmp)) 
    {
      curr.cardinality = sort(table(tmp[j,]),decreasing=TRUE)[[1]];
      curr.score = curr.score + (ncol(tmp) - curr.cardinality)/ncol(tmp);			
    }
    
    score[which(refcol==i)] = 1 - (curr.score/nrow(tmp));
  }
 
 print(score)
 
  
  # Load stage information
  stage = sample(4, ncol(cluster.map), rep=T)
  num.stages = length(unlist(unique(stage)))
  
  # Create annotations
  annotation = data.frame(stage = stage, score=score)
  rownames(annotation) = colnames(cluster.map)
  
  
  # Heatmap settings
  # layout
  fontsize_col = 3
  
  # colors
  # col = diverge_hcl(length(id), c=100) # nice: terrain.colors(5)
  # col= c('#af8dc3', '#f7f7f7', '#7fbf7b')
  # col = terrain.colors(length(id))
  col = brewer.pal(n = length(id), name = 'Set1')
  
  stage.color = brewer.pal(n = num.stages, name = 'YlOrRd')
  score.color = brewer.pal(n = 10, name = 'Greys')
  
  # title
  main= paste("Clustering assignment for subtypes detection (reference ", rownames(cluster.map)[ref.clust],')') 
  
  
  
  
  #####
  annotation_colors = list(stage=stage.color, score=score.color)
  
  
  pheatmap(cluster.map, 
           scale = "none", 
           #clustering_distance_rows = "correlation", 
           cluster_col= F, 
           cluster_rows = F, 
           col=col, 
           main=main,
           fontsize=7,
           fontsize_col=fontsize_col,
           annotation = annotation,
           annotation_colors = annotation_colors,	
           border_color='darkgray',
           border=T,
           margins=c(10,10),
           cellwidth = 3, 
           #cellheight = 18,
           #legend_labels = unlist(unique(id)),
           legend=F,
           #	legend_breaks=-1:4
           #treeheight_column = 20,
  )
  
  return(cluster.map)
  #filename="./results/KO-WT-OV-3way.sig.batch.heatmap.pdf") # name of file output
}

# mapk3 = read.table('vecchi228/nbs_k3.txt')
# mapk4 = read.table('vecchi228/nbs_k4.txt')
# mapk5 = read.table('vecchi228/nbs_k5.txt')

# mapk3 = read.table('oggi228/nbs_k3.txt')
# mapk4 = read.table('oggi228/nbs_k4.txt')
# mapk5 = read.table('oggi228/nbs_k5.txt')

# mapk3 = read.table('oggi220/nbs_k3.txt')
# mapk4 = read.table('oggi220/nbs_k4.txt')
# mapk5 = read.table('oggi220/nbs_k5.txt')

# mapk3 = read.table('oggi208/nbs_k3.txt')
# mapk4 = read.table('oggi208/nbs_k4.txt')
# mapk5 = read.table('oggi208/nbs_k5.txt')

# mapk3 = read.table('fold/nbs_k3.txt')
# mapk4 = read.table('fold/nbs_k4.txt')
# mapk5 = read.table('fold/nbs_k5.txt')

# mapk3 = read.table('1mar/nbs_k3.txt')
# mapk4 = read.table('1mar/nbs_k4.txt')
# mapk5 = read.table('1mar/nbs_k5.txt')

#mapk3 = read.table('ok208/nbs_k3.txt')
#mapk4 = read.table('ok208/nbs_k4.txt')
#mapk5 = read.table('ok208/nbs_k5.txt')
#samples = read.table('ok208/sample_id')

#data = data.frame(mapk3, mapk4, mapk5, row.names = samples[,1])
#colnames(data) = c('K=3', 'K=4', 'K=5')

#cluster.sensitivity(data, reference = 'K=4')
