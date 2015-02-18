#
# oncoPrint : plot a genotype
#



# This is the plotting function
# 
# @param excl.soft A number
# @param col.cluster A number
# @param row.cluster=FALSE
# @param device.new=FALSE 
# @param file=NA
# @param ann.stage=TRUE Show information about stage classification
# @param ann.score=TRUE Show information about the score of a mutation
# @param stage.color='YlOrRd' Color Palette to use with stage
# @param score.color = 'Purples' Color Palette to use with score
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
                      device.new=FALSE, 
                      file=NA, 
                      ann.stage=TRUE, 
                      ann.score=TRUE, 
                      stage.color='YlOrRd', 
                      score.color = 'Purples',  
                      null.color='darkgray', 
                      border.color='white', 
                      font.size=7, 
                      font.column = 3, 
                      font.row = NA, 
                      title= paste('Dataset (genotypes)'),
                      sample.id = FALSE,
                      hide.zeroes = FALSE,
                      legend = TRUE,
                      cellwidth = 5, 
                      cellheigth = 10,
                      group.by.label = FALSE,
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
  
  cat(paste('*** Oncoprint with attributes: stage=', ann.stage, ', score=', ann.score, '\n', sep=''))
  is.compliant(x, 'oncoprint', stage=ann.stage)
  x = enforce.numeric(x)
  
  # We reverse the heatmap under the assumption that ncol(data) << nrow(data)
  data = t(x$genotypes)
  nc = ncol(data) 
  nr = nrow(data)
  
  # If hide.zeros remove event without mutations, that's stronger than a "trim"
  # whihc cuts out the genes with null genotypes
  if (hide.zeroes) {
    cat(paste('Hiding empty genes and samples.\n', sep=''))
    
    data = data[ rowSums(data) != 0, ]
    data = data[ , colSums(data) != 0]
    nr = nrow(data)
    nc = ncol(data)
  }
  
  
  # Sort data, if required. 
  if(excl.sort) {
    cat(paste('Sorting data to enhance exclusivity patterns.\n', sep=''))
    sorted.data = exclusivity.sort(data)
    data = sorted.data$M	
  }	
  
  # If group.by.label group events involving the same label
  if (group.by.label) {
    cat(paste('Grouping events by gene.\n', sep=''))
    genes = as.genes(x)
    data = data[ order(x$annotations[rownames(data), 'event']), ]
  }
  
  cn = colnames(data)
  rn = rownames(data)
  
  # Heatmap score annotation: total 1s per sample
  nmut = colSums(data)
  if(ann.score == TRUE && ann.stage == FALSE) annotation = data.frame(score=nmut)
  if(ann.score == FALSE && ann.stage == TRUE) annotation = data.frame(stage=as.stages(x)[cn, 1])
  if(ann.score == TRUE && ann.stage == TRUE)  annotation = data.frame(stage=as.stages(x)[cn, 1], score=nmut)
  
  if(ann.score == TRUE || ann.stage == TRUE) {
    rownames(annotation) = cn
    annotation_colors = list()
  }
  
  # annotation colors
  if(ann.score == T){
    score.gradient = (colorRampPalette(brewer.pal(6, score.color))) (max(nmut))
    annotation_colors = append(annotation_colors, list(score=score.gradient))
  }
  
  if(ann.stage == T){ 
    different.stages = sort(unique(annotation$stage))
    num.stages = length(different.stages)
    stage.color.attr = append(brewer.pal(n=num.stages, name=stage.color), "#FFFFFF")
    names(stage.color.attr) = append(levels(different.stages), NA)
    annotation_colors = append(annotation_colors, list(stage=stage.color.attr))
  }	  

  # Display also event frequency, which gets computed now
  genes.freq = rowSums(data)/nsamples(x)
  
  # Augment data to make type-dependent colored plots
  types = as.types(x)
  map.gradient = null.color

  for(i in 1:ntypes(x))
  {
    events = as.events(x, type=as.types(x)[i])
    keys = rownames(events)
    sub.data = data[keys, ]

    # shift 1s to 'i', add color to the map  
    idx = which(sub.data == 1)
    if(length(idx) > 0) map.gradient = cbind(map.gradient, as.colors(x)[i])
   
    sub.data[idx] = i
    data[keys, ] = sub.data  
  }
  
  # Augment gene names with frequencies and prepare legend labels
  gene.names = x$annotations[rownames(data),2]
  rownames(data) = paste(round(100 * genes.freq, 0) ,'% ', gene.names, sep='')
  legend.labels = c('0', unique(x$annotations[,1]))
    
  legend.labels = legend.labels[1:(max(data)+1)]
  
  if(is.na(font.row)) 
  {
    font.row = max(c(15 * exp(-0.02 * nrow(data)), 2))    
    cat(paste('Setting automatic font (exp. scaling): ', round(font.row, 1), '\n', sep=''))
  }
  
  # Augment title
  title = paste(title, '\n n = ', nsamples(x),' \t m = ', nevents(x), ' \t |G| = ', ngenes(x),  sep='')
  
  # Pheatmap
  if(ann.score == TRUE || ann.stage == TRUE)  
    pheatmap(data, 
             scale = "none", 
             col = map.gradient, 
             cluster_cols = col.cluster,
             cluster_rows = row.cluster,
             main= title,
             fontsize= font.size,
             fontsize_col= font.column,
             fontsize_row= font.row,
             annotation = annotation,
             annotation_colors = annotation_colors,	
             border_color = border.color,
             border=T,
             margins=c(10,10),
             cellwidth = cellwidth, 
             cellheigth = cellheigth,
             legend=legend,
             legend_breaks = c(0:max(data)),
             legend_labels = legend.labels,
             drop_levels=T,
             show_colnames = sample.id,
             filename=file,
             ...
    )
  else
    pheatmap(data, 
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

