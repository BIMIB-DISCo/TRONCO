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
                      ann.hits = TRUE, 
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
                      cellwidth = NA, 
                      cellheight = NA,
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

  if(excl.sort && nevents(x) > 1) {
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
  	data = data[order(rowSums(data), decreasing = TRUE), ]  	
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
  if(hasGroups) samples.annotation$group = group.samples[cn, 1]

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
	# print(group.color.attr)
	# print(unique(group.samples[,1]))
	# print(samples.annotation)
	
  	names(group.color.attr) = unique(group.samples[,1])
    annotation_colors = append(annotation_colors, list(group=group.color.attr))
   }

    # Augment gene names with frequencies and prepare labels 	
 	 genes.freq = rowSums(data)/nsamples(x)
     gene.names = x$annotations[rownames(data),2]
     gene.names = paste(round(100 * genes.freq, 0) ,'% ', gene.names, sep='') # row labels

	# print(gene.names)

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


		if(length(pathways.color) == 1 && pathways.color %in% rownames(brewer.pal.info))
		{
			cat('Annotating pathways with RColorBrewer color palette', pathways.color, '.\n')
			pathway.colors = append(brewer.pal(n=length(names), name=pathways.color), "#FFFFFF")
		}
		else{
			if(length(pathways.color) != length(names)) 
				stop('You did not provide enough colors to annotate', length(names), 'pathways. 
						Either set pathways.color to a valid RColorBrewer palette or provide the explicit correct number of colors.')
				
			cat('Annotating pathways with custom colors', paste(pathways.color, collapse=','), '.\n')
			pathway.colors = append(pathways.color, "#FFFFFF")
		}
		names(pathway.colors) = append(names, NA)

		pathway.colors = pathway.colors[ unique(genes.annotation$pathway) ]
		
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
      sub.data = data[keys.subset, , drop = FALSE]
      
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
             # margins=c(10,10),
             cellwidth = cellwidth, 
             cellheight = cellheight,
             legend = legend,             
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
    
    return(ret)
}


##### Pathway print
pathway.visualization = function(x, 
	title = paste('Pathways:', paste(names(pathways), collapse=', ', sep='')), 
	file, 
	pathways.color = 'Set2', 
	aggregate.pathways, 
	pathways,
	...) 
{	
	names = names(pathways)
	    
    if(length(pathways.color) == 1 && pathways.color %in% rownames(brewer.pal.info))
	{
		cat('Annotating pathways with RColorBrewer color palette', pathways.color, '.\n')
		pathway.colors = brewer.pal(n=length(names), name=pathways.color)
	}
	else{
		print(pathways.color)		
		if(length(pathways.color) != length(names)) 
			stop('You did not provide enough colors to annotate ', length(names), ' pathways. 
					Either set pathways.color to a valid RColorBrewer palette or provide the explicit correct number of colors.')
				
		cat('Annotating pathways with custom colors', paste(pathways.color, collapse=','), '.\n')
		pathway.colors = pathways.color
	}
	
	names(pathway.colors) = names
	
    cat(paste('*** Processing pathways: ', paste(names, collapse=', ', sep=''), '\n', sep=''))
  
	cat(paste('\n[PATHWAY \"', names[1],'\"] ', paste(pathways[[1]], collapse=', ', sep=''), '\n', sep=''))

    data.pathways = as.pathway(x, pathway.genes=pathways[[1]], 
                             pathway.name=names[1], aggregate.pathway = aggregate.pathways)
	
	data.pathways = change.color(data.pathways, 'Pathway', pathways.color[1])
	data.pathways = rename.type(data.pathways, 'Pathway', names[1])

  if(length(names) > 1)
  {  
    for(i in 2:length(pathways))
  	{
  	  cat(paste('\n\n[PATHWAY \"', names[i],'\"] ', paste(pathways[[i]], collapse=', ', sep=''), '\n', sep=''))
  	  pathway = as.pathway(x, pathway.genes=pathways[[i]], pathway.name=names[i], aggregate.pathway = aggregate.pathways)
                                       
	  pathway = change.color(pathway, 'Pathway', pathway.colors[i])
      pathway = rename.type(pathway, 'Pathway', names[i])
     
      # show(pathway)
      # # print(has.stages(pathway))
       # print('bindo..')
	  # print(has.stages(data.pathways))

  	  data.pathways = ebind(data.pathways, pathway)
      # show(data.pathways)
  	}
  }

  # data.pathways = enforce.numeric(data.pathways)
 # show(data.pathways)
  ret = oncoprint(trim(data.pathways), title=title, file=file, ...)

  return(ret)
}


# #### Consensus matrix (intra-clusters)
# oncoprint.consensus = function(models, MIN.HITS = 0)
# {


# smaller.to.bigger = function(m,cn)
# {
	# x = matrix(0, nrow = length(cn), ncol = length(cn))
	# rownames(x) = cn
	# colnames(x) = cn
 
	
	# for(i in 1:nrow(m))
		# for(j in 1:nrow(m))
			# x[rownames(m)[i], rownames(m)[j]] = ifelse(m[i,j] == 1, 1, 0) 
	# return(x)
# }


	# # All the adjacency matrices
	# matrices = list()
	# for(i in 1:length(models))
		# matrices = append(matrices, list(models[[i]]$adj.matrix$adj.matrix.bic))

	# # All their colnames - all possible events and types
	# cn = unique(Reduce(union, lapply(matrices, colnames)))
	
	# all.events = NULL
	# for(i in 1:length(models)) all.events = rbind(all.events, as.events(models[[i]]$data))
	# all.events = unique(all.events)

	# all.types = NULL
	# for(i in 1:length(models)) all.types = rbind(all.types, models[[i]]$data$types)
	# all.types = unique(all.types)
	
	# # Consensus + overall adjacency matrix
	# consensus = Reduce('+', lapply(matrices, smaller.to.bigger, cn=cn))
	# adjacency = consensus
	# adjacency[adjacency < MIN.HITS] = 0 
	# adjacency[adjacency > 1] = 1 
		
	# cat('Consensus adjacency matrix:', nrow(adjacency), 'x', ncol(adjacency), ', minimum consensus', MIN.HITS, '\n')
	
	
			
	# keys = Reduce(rbind, lapply(models, as.events))
	
	# types = data.frame(event=rep('pattern', nrow(super.model)), 
		# stringsAsFactors = F)
	
	# labels = data.frame(name=rep('none', nrow(super.model)), 
		# stringsAsFactors = F)
		
	# pattern.type =	data.frame(pattern=rep('NA', nrow(super.model)), 
		# stringsAsFactors = F)


	# for(i in 1:length(rownames(super.model)))
	# {
		# types$event[i] = keys[rownames(super.model)[i], 'type']		
				
		# labels$name[i] = keys[rownames(super.model)[i], 'event']				
	
		# prefix = gsub("_.*$", "", labels$name[i])
		# prefix = gsub( "\"","", prefix)
		
		# if(prefix %in% c('AND', 'OR', 'XOR'))
		# {			
			# pattern.type$pattern[i] = prefix
			# if(prefix == 'AND') pattern.type$pattern[i] = 'co-occurrence'
			# if(prefix == 'OR') pattern.type$pattern[i] = 'soft exclusivity'
			# if(prefix == 'XOR') pattern.type$pattern[i] = 'hard exclusivity'
		
			# compact.label = strsplit(labels$name[i], '_' )[[1]]
			# labels$name[i] = paste(
				# compact.label[2:length(compact.label)],
				# collapse = ' / ')
		# }
		# rownames(super.model)[i] = 
			# paste(
				# keys[rownames(super.model)[i], 'event'], 
				# keys[rownames(super.model)[i], 'type'])  
	# }
	
	# colnames(super.model) = rownames(super.model)
	# rownames(types) = rownames(super.model)
	# rownames(labels) = rownames(super.model)
	# rownames(pattern.type) = rownames(super.model)

	# types[types == 'Hypothesis'] = 'Pattern'

	# ########### Filetering for recurrrent
	# view = super.model
	# view = view[apply(view, 1, max) > MIN,]
	# view = view[, apply(view, 2, max) > MIN]
	
	# ########### Ordering according to the number of times a relation is found
	# group.order = function(x, k, margin)
	# {		
		# ord = which(apply(x, margin, max) == k)
		
		# if(margin == 1)
		# {
			# to.sort = x[ord, , drop = F]
			# to.sort = to.sort[order(rowSums(to.sort), decreasing = T), , drop = F]
			# rownames(x)[ord] = rownames(to.sort)		
			# x[ord, ] = to.sort
		# }
		# else
		# {
			# to.sort = x[, ord , drop = F]
			# to.sort = to.sort[, order(colSums(to.sort), decreasing = T), drop = F]
			# colnames(x)[ord] = colnames(to.sort)
			# x[, ord] = to.sort	
		# }
		
		# return(x)
	# }	
		
	# selects.recurr = order(apply(view, 1, max), decreasing = T)
	# view = view[selects.recurr, ]	
	
	# max.values = apply(view, 1, max)
	# for(i in 1:length(unique(max.values)))
	 	# view = group.order(view, unique(max.values)[i], 1)
	
	# selected.recurr = order(apply(view, 2, max), decreasing = T)
	# view = view[, selected.recurr]	

	# max.values = apply(view, 2, max)
	# for(i in 1:length(unique(max.values)))
	 	# view = group.order(view, unique(max.values)[i], 2)

	# gaps.row = match(
		# unique(apply(view, 1, max)),
		# apply(view, 1, max)) - 1

	# gaps.col = match(
	 	# unique(apply(view, 2, max)),
		# apply(view, 2, max)) - 1


	# ############
	# view.ones = view
	# view.ones[view.ones > 1] = 1
	# ann.row = data.frame(selects=rowSums(view.ones), 
		# row.names = rownames(view),
		# stringsAsFactors = F)
	
	# # ann.row$selects = .bincode(ann.row$selects,
		# # c(0, 
		# # max(ann.row$selects) / 4,
		# # max(ann.row$selects) / 2,
		# # max(ann.row$selects) * 3 / 4,
		# # max(ann.row$selects)), 
		# # TRUE)		

	
	# ann.row = cbind(
		# types[rownames(ann.row), , drop = FALSE],	
		# pattern.type[rownames(ann.row), , drop = FALSE],
		# ann.row)

	# ann.col = data.frame(selected=colSums(view.ones), 
		# row.names = colnames(view),
		# stringsAsFactors = F)

	# ann.col = cbind(
		# types[rownames(ann.col), , drop = FALSE],	
		# pattern.type[rownames(ann.col), , drop = FALSE],
		# ann.col)
	
	# selection = max(max(ann.row$selects), max(ann.col$selected))
	# selection.palette = brewer.pal(9, 'PuBuGn')[4:9]
	# # selection.palette = colorRampPalette(selection)(selection)
	
	# ann.colors = list(
		# event = as.colors(Cluster3.methylation_subtype.aic$data),
		# pattern = c('red', 'darkgreen', 'orange', 'white'),
		# # selects = colorRampPalette(brewer.pal(9, 'PuBuGn'))(max(ann.row$selects) + 1),
		# # selected = colorRampPalette(brewer.pal(9, 'PuBuGn'))(max(ann.col$selected) + 1)
		# # selects = selection.palette,
		# # selected = selection.palette
		# selects = colorRampPalette(selection.palette)(max(ann.row$selects) + 1),
		# selected = colorRampPalette(selection.palette)(max(ann.col$selected) + 1)
		# )	
	# ann.colors$event['Hypothesis'] = 'white'
	# names(ann.colors$event)[4] = 'Pattern'
	# names(ann.colors$pattern) = c('hard exclusivity', 'co-occurrence', 'soft exclusivity', 'NA')
	
	# print(ann.colors$selects)
		
	# map.color = colorRampPalette(brewer.pal(3, 'YlOrBr'))(max(super.model) + 1)
	# # map.color = brewer.pal(max(super.model) + 1, 'YlOrBr')
	# map.color[1] = 'gray92'
	
	# map.nozeroes = view
	# map.nozeroes[view == 0] = ''
	
	# pheatmap(view, 
		# main = 'Selectivity relations (intra-clusters consensus)',
		# scale = 'none', 
		# color = map.color,
		# font.size=4,
		# fontsize_col = 6,
		# fontsize_row = 6,
		# annotation_row = ann.row,
		# annotation_col = ann.col,
		# annotation_colors = ann.colors,
		# fontsize_number = 6,
		# cluster_rows = F,
		# cluster_cols = F,
		# labels_row = labels[rownames(view), ],
		# labels_col = labels[colnames(view), ],
		# number_format = '%d',
		# display_numbers = map.nozeroes,
		# gaps_col = gaps.col,
		# gaps_row = gaps.row,
		# border_color = 'lightgray'
		# )
# }	


