#' Return all events involving certain genes and types
pattern.events = function(x, hypothesis)
{
	hevents = x$hypotheses$patterns[[hypothesis]]
	if(length(hevents)==0){
		stop('The hypothesis is not valid!')
	}
	return(hevents)
}

# Return the names of the patterns in the dataset
as.patterns = function(x)
{
	return(names(x$hypotheses$patterns))
}



pattern.plot = function(x, hypo, to)
{
	# keys
	keys = pattern.events(x, hypo)
	events = as.events(x)
	
	cat('Pattern elements.\n')
	events.names = events[keys, , drop = FALSE]
	print(events.names)
	

	cat('Pattern tested against:', to[1], to[2], '\n')
		

	# print(keys)
	# print(events)
	
	
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
	
	
	# 
	# print(to.samples)
	
	
	# print(mat)
	
	# Choose colors according to event type
	grid.col = rep('gray', nrow(matrix) + 1) 
	names(grid.col) = c(rownames(matrix), colnames(matrix))
	for(i in 1:length(keys))
		grid.col[keys[i]] = as.colors(x)[events.names[keys[i], 'type' ]]
	grid.col[colnames(matrix)] = as.colors(x)[as.events(x, genes = to[1], types=to[2])[, 'type' ]]
	grid.col['soft'] = 'orange'
		
	
	# Informative labels
	for(i in 1:length(keys))
	{
		# spaces = nevents(x, genes= events.names[i, 'event' ])
		
		rownames(matrix)[i] = paste(events.names[i, 'event' ], paste(rep(' ', i), collapse=''))
		names(grid.col)[i] = rownames(matrix)[i]
		
		# rownames(matrix)[i] = paste(events.names[i, 'event' ], events.names[i, 'type' ])
		# names(grid.col)[i] = paste(events.names[i, 'event' ], events.names[i, 'type' ])
	}

	print(grid.col)
	print(matrix)

	# gap.degree = c(rep(2, 4), 10, rep(2, 4), 10) 
	# circos.clear() 
	# circos.par(gap.degree = gap.degree, start.degree = -10/2
	circos.clear()
	circos.par(gap.degree=rep(c(2, 10), 7))
	
	chordDiagram(matrix, 
		grid.col = grid.col,
		annotationTrack = "grid", 
		preAllocateTracks = list(track.height = 0.3)
		)
	
	circos.trackPlotRegion(
		track.index = 1, 
		panel.fun = function(x, y) 
			{ 
				xlim = get.cell.meta.data("xlim") 
				ylim = get.cell.meta.data("ylim") 
				sector.name = get.cell.meta.data("sector.index") 
				circos.text(mean(xlim), ylim[1], sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5)) 
				},
		 bg.border = NA)	
	
}
















