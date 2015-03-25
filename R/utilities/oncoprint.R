# oncoPrint : plot a genotype
# 
# @param excl.soft A number
# @param col.cluster A number
# @param row.cluster=FALSE
# @param device.new=FALSE 
# @param file=NA
# @param ann.stage=TRUE Show information about stage classification
# @param ann.hits=TRUE Show information about the hits in a gene
# @param stage.color='YlOrRd' Color Palette to use with stage
# @param hits.color = 'Purples' Color Palette to use with score
# @param null.color='darkgray' Background color
# @param border.color='white' 
# @param font.size=7
# @param font.column = 3 
# @param title= paste('Genotypes')
# @param sample.id = F Show sample name at the bottom of the heatmap
# @param hide.zeroes = F Hide events without mutations
oncoprint <- function(x, 
                      excl.sort=TRUE, 
                      col.cluster=FALSE, 
                      row.cluster=FALSE, 
                      file=NA, 
                      ann.stage = has.stages(x), 
                      ann.hits=TRUE, 
                      stage.color='YlOrRd', 
                      hits.color = 'Purples',  
                      null.color='lightgray', 
                      border.color='white', 
                      font.size=7, 
                      font.column = 3, 
                      font.row = NA, 
                      title= paste('Dataset (genotypes)'),
                      sample.id = FALSE,
                      hide.zeroes = FALSE,
                      legend = TRUE,
                      legend.cex = 1.0,
                      cellwidth = 5, 
                      cellheigth = 10,
                      group.by.label = FALSE,
                      group.samples = NA,
                      pathways = NA,
                      pathways.color = 'Set1',
                      ...) 
{
  if (!require('pheatmap')) {
    install.packages('pheatmap', dependencies = TRUE)
    library(pheatmap)
  }
  
  if (!require('RColorBrewer')) {
    install.packages('RColorBrewer', dependencies = TRUE)
    library(RColorBrewer)
  }
  
  # This function sorts the matrix for better visualization of mutual exclusivity across genes
  exclusivity.sort <- function(M) {
    geneOrder <- sort(rowSums(M), decreasing=TRUE, index.return=TRUE)$ix;
    scoreCol <- function(x) {
      score <- 0;
      for(i in 1:length(x)) {
        if(x[i]) {
          score <- score + 2^(length(x)-i);
        }
      }
      return(score);
    }
    scores <- apply(M[geneOrder, ], 2, scoreCol);
    sampleOrder <- sort(scores, decreasing=TRUE, index.return=TRUE)$ix;
    
    res = list()
    res$geneOrder = geneOrder
    res$sampleOrder = sampleOrder
    res$M = M[geneOrder, sampleOrder]
    
    return(res);
  }
  
  cat(paste('*** Oncoprint with attributes: stage=', ann.stage, ', hits=', ann.hits, '\n', sep=''))
  is.compliant(x, 'oncoprint', stage=ann.stage)
  x = enforce.numeric(x)
  
  # If hide.zeros trim x
  if (hide.zeroes) {
    cat(paste('Trimming the input dataset (hide.zeroes).\n', sep=''))
   	x = trim(x)
   	}
 
  # We reverse the heatmap under the assumption that ncol(data) << nrow(data)
  data = t(x$genotypes)
  nc = ncol(data) 
  nr = nrow(data)
    
  # Sort data, if required. excl.sort and group.samples are not compatible
  hasGroups = !any(is.na(group.samples))

  if(excl.sort && hasGroups) 
  	stop('Disable sorting for mutual exclusivity (excl.sort=FALSE) or avoid using grouped samples (group.samples=NA).')

  if(excl.sort) {
    cat(paste('Sorting samples ordering to enhance exclusivity patterns.\n', sep=''))
    sorted.data = exclusivity.sort(data)
    data = sorted.data$M	
  }
  
  if(hasGroups)
  {
  	grn = rownames(group.samples)
  	
    cat(paste('Grouping samples according to input groups (group.samples).\n', sep=''))
    
    if(any(is.null(grn))) stop('Input groups should have sample names.')
    
  	if(!setequal(grn, as.samples(x)))
  		stop(paste0('Missing group assignment for samples: ', paste(setdiff(as.samples(x), grn), collapse=', '),'.'))
 
 	# Order groups by label, and then data (by column)
 	order = order(group.samples)
  	group.samples = group.samples[order, , drop=FALSE]
  	
  	data = data[, rownames(group.samples)]  	  	
  }	
  
  ##### If group.by.label group events involving the gene symbol
  if (group.by.label) {
    cat(paste('Grouping events by gene.\n', sep=''))
    genes = as.genes(x)
    data = data[ order(x$annotations[rownames(data), 'event']), ]
  }
  
  cn = colnames(data)
  rn = rownames(data)
  
  ##### SAMPLES annotations: hits (total 1s per sample), stage or groups
  samples.annotation = NA
  nmut = colSums(data)
  
  if(ann.hits == TRUE && ann.stage == FALSE) samples.annotation = data.frame(hits=nmut)
  if(ann.hits == FALSE && ann.stage == TRUE) samples.annotation = data.frame(stage=as.stages(x)[cn, 1])
  if(ann.hits == TRUE && ann.stage == TRUE)  samples.annotation = data.frame(stage=as.stages(x)[cn, 1], hits=nmut)
  if(hasGroups) samples.annotation $group = as.factor(group.samples[cn, 1])

  ##### Color each annotation 
  if(ann.hits || ann.stage || hasGroups) {
    rownames(samples.annotation) = cn
    annotation_colors = list()
  }

  if(ann.hits){
    hits.gradient = (colorRampPalette(brewer.pal(6, hits.color))) (max(nmut))
    annotation_colors = append(annotation_colors, list(hits=hits.gradient))
  }
  
  if(ann.stage){ 
    different.stages = sort(unique(samples.annotation$stage))
    num.stages = length(different.stages)
    stage.color.attr = append(brewer.pal(n=num.stages, name=stage.color), "#FFFFFF")
    names(stage.color.attr) = append(levels(different.stages), NA)
    annotation_colors = append(annotation_colors, list(stage=stage.color.attr))
  }
  
  if(hasGroups)	{
  	ngroups = length(unique(group.samples[,1]))
  	group.color.attr = brewer.pal(n=ngroups, name='Accent')
  	names(group.color.attr) = levels(as.factor(unlist(unique(group.samples[,1]))))
    annotation_colors = append(annotation_colors, list(group=group.color.attr))
   }
   
    # Augment gene names with frequencies and prepare labels 	
 	 genes.freq = rowSums(data)/nsamples(x)
     gene.names = x$annotations[rownames(data),2]
     gene.names = paste(round(100 * genes.freq, 0) ,'% ', gene.names, sep='') # row labels


	# GENES ANNOTATIONS - PATHWAYS
	genes.annotation = NA

    if(!all(is.na(pathways)))
    {
		names = names(pathways)  	
				
		genes.annotation = data.frame(row.names = rn, stringsAsFactors = FALSE)
		genes.annotation$pathway = rep(NA, nrow(data))
		
		for(i in 1:length(names)) 
		{
			pathway = names[i]
			genes.pathway = rownames(as.events(x, genes=pathways[[names[i]]]))
			genes.annotation[genes.pathway, 'pathway'] = names[i] 
		}

		# print(annotation_colors)			

		pathway.colors = append(brewer.pal(n=length(names), name=pathways.color), "#FFFFFF")
		names(pathway.colors) = append(names, NA)

		annotation_colors = append(annotation_colors, list(pathway=pathway.colors))
		# print(annotation_colors)				   	
   }   
  
  # Augment data to make type-dependent colored plots
  types = as.types(x)
  map.gradient = null.color
  
  for(i in 1:ntypes(x))
  {
    events = as.events(x, type=as.types(x)[i])
    keys = rownames(events)
    
    if (ntypes(x) > 1) {
      keys.subset = keys[unlist(lapply(keys, function(x, data){if (x %in% data) T else F}, rownames(data)))]
      sub.data = data[keys.subset, ]
      
      # shift 1s to 'i', add color to the map  
      idx = which(sub.data == 1)
      if(length(idx) > 0) map.gradient = cbind(map.gradient, as.colors(x)[i])
     
      sub.data[idx] = i
      data[keys.subset, ] = sub.data 

    } else {
      map.gradient = cbind(map.gradient, as.colors(x)[i])
    }
  }
    
  if(is.na(font.row)) 
  {
    font.row = max(c(15 * exp(-0.02 * nrow(data)), 2))    
    cat(paste('Setting automatic font (exponential scaling): ', round(font.row, 1), '\n', sep=''))
  }
  
  # Augment title
  title = paste(title, '\n n = ', nsamples(x),'    m = ', nevents(x), '    |G| = ', ngenes(x),  sep='')
  
    legend.labels = c('none', unique(x$annotations[,1]))
    
     legend.labels = legend.labels[1:(max(data)+1)]

  
  # Pheatmap
  if(ann.hits == TRUE || ann.stage == TRUE || hasGroups)  
   ret = pheatmap(data, 
             scale = "none", 
             col = map.gradient, 
             cluster_cols = col.cluster,
             cluster_rows = row.cluster,
             main= title,
             fontsize = font.size,
             fontsize_col= font.column,
             fontsize_row= font.row,
             annotation_col = samples.annotation,
             annotation_row = genes.annotation,
             annotation_colors = annotation_colors,	
             border_color = border.color,
             border=T,
             #margins=c(10,10),
             cellwidth = cellwidth, 
             cellheigth = cellheigth,
             legend=legend,
             legend_breaks = c(0:max(data)),
             legend_labels = legend.labels,
             legend.cex = legend.cex,
             labels_row = gene.names,
             drop_levels=T,
             show_colnames = sample.id,
             filename=file,
             ...
    )
  else
    ret = pheatmap(data, 
             scale = "none", 
             col = map.gradient, 
             cluster_cols = col.cluster,
             cluster_rows = row.cluster,
             main= title,
             fontsize= font.size,
             fontsize_col= font.column,
             fontsize_row= font.row,
             border_color= border.color,
             border=T,
             margins=c(10,10),
             cellwidth = cellwidth, 
             cellheigth = cellheigth,
             legend=legend,
             legend_breaks = c(0:max(data)),
             legend_labels = legend.labels,
             show_colnames = sample.id,
             filename=file,
             ...
    )
    
}


##### Pathway print
pathway.visualization = function(x, file, aggregate.pathways, names, ...) 
{	
	input = list(...)

  if(length(names) != length(input))
    stop('Missing pathway names...')
  
	#name=deparse(substitute(data)),
  cat(paste('*** Processing pathways: ', paste(names, collapse=', ', sep=''), '\n', sep=''))
  
	cat(paste('\n[PATHWAY \"', names[1],'\"] ', paste(unlist(input[1]), collapse=', ', sep=''), '\n', sep=''))
  data.pathways = as.pathway(x, pathway.genes=unlist(input[1]), 
                             pathway.name=names[1], aggregate.pathway = aggregate.pathways)
	
  if(length(names) > 1)
  {  
    for(i in 2:length(input))
  	{
  	  cat(paste('\n\n[PATHWAY \"', names[i],'\"] ', paste(unlist(input[i]), collapse=', ', sep=''), '\n', sep=''))
  	  data.pathways = ebind(data.pathways,
                            as.pathway(x, pathway.genes=unlist(input[i]), 
                                       pathway.name=names[i], aggregate.pathway = aggregate.pathways))
  	}
  }


  oncoprint(trim(data.pathways), title=paste('Pathways:', paste(names, collapse=', ', sep='')), 
          file=file)

  return(data.pathways)
}

