as.samples = function(x)
{
	return(rownames(x$genotypes))
}

as.genes = function(x)
{
	return(unlist(unique(x$annotations[, 'event'])))
}

as.gene = function(x, g)
{
	data = data.frame(x$genotypes[, which(as.events(x)[, 'event'] == g)], row.names = as.samples(x))
	colnames(data) = x$annotations[which(as.events(x)[, 'event'] == g)]
	
	return(data)
}

as.events = function(x)
{
	return(x$annotations[,c('type', 'event')])
}

as.stages = function(x)
{
	if(!is.null(x$stages))
	return(x$stages)
}

as.types = function(x)
{
	return(rownames(x$types))
}

as.colors = function(x)
{
	return(x$types[, 'color'])
}

show = function(x, view = 10)
{
	is.compliant(x)
	
	cat(paste('Dataset: n=', nsamples(x), ', m=', nevents(x), ', |G|=', ngenes(x), '.\n', sep=''))
	cat(paste('Events: ', paste(as.types(x), collapse=', '), '.\n', sep=''))
	cat(paste('Colors: ', paste(as.colors(x), collapse=', '), '.\n', sep=''))

	if(!is.null(x$stages))
	{
		cat(paste('Stages: '))
		# cat(paste(paste('\t', rownames(as.stages(x)[1: view,]), ':', as.stages(x)[1: view, 1], sep=' '), collapse='\n'))
		
		s = unlist(unique(as.stages(x)[, 1]))
		cat(paste(s, collpase=', '), '.\n', sep='')
	 }

	cat(paste('Events: ', view, ' shown.\n'))
	cat(paste(paste('\t', rownames(as.events(x)[1: view,]), ':', as.events(x)[1: view, 1], as.events(x)[1: view, 2], sep=' '), collapse='\n'))

	cat(paste('\nGenotypes for ', view, ' events.\n'))
	head(x$genotypes[,1:view])
	
	# head(x$annotations[1:nev,])

}


nsamples = function(x)
{
	return(nrow(x$genotypes))
}

nevents = function(x)
{
	return(ncol(x$genotypes))
}

ngenes = function(x)
{
	return(length(as.genes(x)))
}



