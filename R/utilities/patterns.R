#' Return all events involving certain genes and types
pattern.events = function(x, hypothesis)
{
	#	hevents = x$patterns[[hypothesis]]
	if(is.null(x$patterns)) 
		hevents = x$hypotheses$patterns[[hypothesis]]
	else
		hevents = x$patterns[[hypothesis]]

	if(length(hevents)==0){
		stop('The hypothesis is not valid!')
	}
	return(hevents)
}

# Return the names of the patterns in the dataset
as.patterns = function(x)
{
	# return(names(x$hypotheses$patterns))
	if(is.null(x$patterns)) return(names(x$hypotheses$patterns))

	return(names(x$patterns))
}




## SPOSTARE IN UN FILE APPOSITO?
#
#
pattern.plot = function(x, group, to, gap.cex = 1.0, legend.cex = 1.0, label.cex = 1.0, title=paste(to[1], to[2] ))
{
	# keys
	events = as.events(x)
	
	keys = NULL
	for(i in 1:nrow(group))
		keys = c(keys, rownames(as.events(x, genes=group[i, 'event'], types=group[i, 'type'])))

	cat('Group:\n')
	events.names = events[keys, , drop = FALSE]
	print(events.names)
	

	cat('Group tested against:', to[1], to[2], '\n')
	
	library(circlize)
	
	# HARD exclusivity: 1 for each pattern element
	# CO-OCCURRENCE: 1
	# SOFT exclusivity: 1
	# OTHERS: 1 
	matrix = matrix(0, nrow = length(keys) + 3, ncol = 1)
	rownames(matrix) = c(keys, 'soft', 'co-occ', 'other')
	# colnames(matrix) = paste(to, collapse=':')
	colnames(matrix) = to[1]
	
	to.samples = which.samples(x, gene=to[1], type=to[2])
	cat('Pattern conditioned to samples:', paste(to.samples, collapse=', '), '\n')
	
	pattern.samples = list()
	hard.pattern.samples = list()
	for(i in 1:length(keys))
	{
		# Samples
		samples = which.samples(x, gene= events.names[i, 'event'], type= events.names[i, 'type'])
		
		# Pattern samples
		pattern.samples = append(pattern.samples, list(samples))
		
		# Other negative samples
		negative.samples = list()
		for(j in 1:length(keys))
			if(i != j)
				negative.samples = append(negative.samples, 
					list(which.samples(x, gene= events.names[j, 'event'], type= events.names[j, 'type'], neg = TRUE)))
		
		pattern.negative = Reduce(intersect, negative.samples)
		hard.pattern.samples = append(hard.pattern.samples, list(intersect(samples, pattern.negative)))		
	}
	names(pattern.samples) = keys
	names(hard.pattern.samples) = keys
	
	# print(hard.pattern.samples)
	
	# print('tutti:')
	# print(unlist(pattern.samples))

	# CO-OCCURRENCES	
	pattern.co.occurrences = Reduce(intersect, pattern.samples)
	co.occurrences = intersect(to.samples, pattern.co.occurrences)
	matrix['co-occ', ] = length(co.occurrences)
	cat('Co-occurrence in #samples: ', matrix['co-occ', ], '\n')

	# HARD EXCLUSIVITY		
	for(i in 1:length(keys))
		matrix[keys[i], ] = length(intersect(to.samples, hard.pattern.samples[[keys[i]]])) 
	cat('Hard exclusivity in #samples:', matrix[keys, ], '\n')	
	
	# OTHER
	union = unique(unlist(pattern.samples))
	union = setdiff(as.samples(x), union)
	matrix['other', ] = length(intersect(to.samples, union))
	cat('Other observations in #samples:', matrix['other', ], '\n')	

	# SOFT EXCLUSIVITY
	matrix['soft', ] = length(to.samples) - colSums(matrix)
	cat('Soft exclusivity in #samples:', matrix['soft', ], '\n')	
	
	# Choose colors according to event type
	sector.color = rep('gray', nrow(matrix) + 1) 
	link.color = rep('gray', nrow(matrix)) 
	
	names(sector.color) = c(rownames(matrix), colnames(matrix))
	names(link.color) = rownames(matrix)

	add.alpha <- function(col, alpha=1)
	{
  	if(missing(col)) stop("Please provide a vector of colours.")
  	apply(sapply(col, col2rgb)/255, 2, function(x) rgb(x[1], x[2], x[3], alpha=alpha))  
     }

	# Link colors - event types	
	for(i in 1:length(keys))
		link.color[keys[i]] = as.colors(x)[events.names[keys[i], 'type' ]]
		
	# print(link.color)
		
	link.color['soft'] = 'orange'
	link.color['co-occ'] = 'darkgreen'
	
	
	idx.max = which(matrix == max(matrix))
	link.style = matrix(0, nrow=nrow(matrix), ncol=ncol(matrix))
	rownames(link.style) = rownames(matrix)
	colnames(link.style) = colnames(matrix)
	link.style[idx.max, 1] = 5
	print(link.style)
	
	# link.color[idx.max] = add.alpha(link.color[idx.max], .3)
	# print(link.color)
	
		
	# Sector colors - event types	
	sector.color[1:length(keys)] = 'red' # XOR
	sector.color['soft'] = 'orange'      # OR
	sector.color['co-occ'] = 'darkgreen' # AND
	sector.color['co-occ'] = 'darkgreen' # targer
	sector.color[colnames(matrix)] = as.colors(x)[as.events(x, genes = to[1], types=to[2])[, 'type' ]]
	
	# print(link.color)
	# print(sector.color)
		
	# Informative labels
	for(i in 1:length(keys))
	{
		# rownames(matrix)[i] = paste(paste(rep(' ', i), collapse=''), events.names[i, 'event' ])
		if(nevents(x, genes = events.names[i, 'event' ]) > 1)
			rownames(matrix)[i] = paste(paste(rep(' ', i), collapse=''), events.names[i, 'event' ])
		else rownames(matrix)[i] = events.names[i, 'event' ]
		
		names(sector.color)[i] = rownames(matrix)[i]		
	}

	cat('Circlize matrix.\n')
	print(matrix)

	circos.clear()
	
	gaps = c(rep(2, length(keys) - 2), rep(15 * gap.cex, 4), rep(40 * gap.cex, 2))
	
	circos.par(gap.degree = gaps)

	# matrix = matrix[ order(matrix[, 1]) , , drop = FALSE]	
		# matrix = matrix[ order(matrix[, 1], decreasing = T) , , drop = FALSE]	
		# print(matrix)
		
	
		
		
	chordDiagram(matrix, 
		grid.col = sector.color,
		annotationTrack = "grid", 
		preAllocateTracks = list(track.height = 0.3),
		row.col = link.color,
		link.border = 'black',
		link.lty = link.style,
		link.lwd = 0.3
		)
		
	for(si in get.all.sector.index()) 
	{ 
		# here the index for the grid track is 2 
		circos.axis(h = "top", labels.cex = 0.3, major.tick.percentage = .4, sector.index = si, track.index = 2) 
	}	
	
	circos.trackPlotRegion(
		track.index = 1, 
		panel.fun = function(x, y) 
			{ 
				xlim = get.cell.meta.data("xlim") 
				ylim = get.cell.meta.data("ylim") 
				sector.name = get.cell.meta.data("sector.index") 
				circos.text(mean(xlim), cex = 1.0 * label.cex, ylim[1], sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5)) 
				},
		 bg.border = NA)	
	
	group.legend = c(
		expression(bold('Group')),
		as.vector(unique(events.names[, 'event'])),
		expression(bold('Target')),
		to[1]
	)
	
	print(group.legend)
	
	legend('bottomleft', 
		cex = .6 * legend.cex, 
		pt.cex = 1, 
		bty='n',
		legend = group.legend
		)
		
	events.legend.pch = rep(19, length(unique(events.names[, 'type'])))	
		
	legend('bottomright',
           legend = unique(events.names[, 'type']),
           title = expression(bold('Events type')),
           bty = 'n',
           cex = .6 * legend.cex,
           pt.cex = 1,
           pch = events.legend.pch,
           col = as.colors(x)[unique(events.names[, 'type'])]
           # pt.bg = pt_bg
           )
		
	legend("topright", 
		cex = 0.6  * legend.cex, 
		pt.cex = 1, 
		title = expression(bold('Observations')), 
		horiz=F, 
		bty='n', 
		legend = c(
				paste(sum(matrix[1:length(keys) ,]), 'hard exclusivity'), 
				paste(matrix[length(keys) + 1, ],  'soft exclusivity'),
				paste(matrix[length(keys) + 2, ], 'co-occurrence'),
				paste(matrix[length(keys) + 3, ], 'none')
				),
			fill= c('red', 'orange', 'darkgreen', 'gray')
			)
	
	if(!is.na(title)) title(title, cex.main = 1.0 * label.cex)
}
















