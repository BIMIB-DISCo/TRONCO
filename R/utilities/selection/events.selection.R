# events.selection - select a subset of the input genotypes 'x'. Selection can be done by frequency and gene symbols.
#
# @x: input genotypes
# @filter.freq: [0,1] value which constriants the minimum frequence of selected events
# @filter.in.names: gene symbols which will be included
# @filter.out.names: gene symbols which will NOT be included
"events.selection" <- function(x, filter.freq=NA, filter.in.names=NA, filter.out.names=NA) {

	is.compliant(x, err.fun='events.selection: input')
	dataset = x$genotypes
	
	cat(paste('*** Events selection: #events=', ncol(dataset), ', #types=', nrow(x$types),  sep=''))
	cat(paste('\nFrequency Filter : ', !is.na(filter.freq), '\n', sep=''))
	cat(paste('In-Names Filter : ', !any(is.na(filter.in.names)), '\n', sep=''))
	cat(paste('Out-Names Filter : ', !any(is.na(filter.out.names)), sep=''))

	if(is.na(filter.out.names) && is.na(filter.in.names) && is.na(filter.freq))
		return(x);

	valid = rep(FALSE, ncol(x$genotypes))
					
	if(!is.na(filter.freq))
	{	
		cat(paste('\nGenes with a minimum mutation frequency of ', filter.freq, ' (', 
              round(nrow(x$genotypes) * filter.freq, 0),' hits) are: ', sep=''))
    enforce.numeric(x)		
    
		for(i in 1:nevents(x))
		{		
			mut.freq <- sum(x$genotypes[,i])/nsamples(x)
			valid[i] <- mut.freq > filter.freq
					
			if(valid[i]) 
			{
				name = x$annotations[colnames(x$genotypes)[i], ]
				cat(paste('\"', name[2], '\" [', name[1],' ', round(mut.freq, 1), '] \n', sep=''))
			}
		}
	}
						
	if(!any(is.na(filter.in.names)))
	{
		cat(paste('\nEvents for the following genes will be selected (filter.in.names): ', 
              sep='', paste(filter.in.names, collapse=', ')))
		
		colnames = which(x$annotations[,2] %in% filter.in.names, arr.ind=T)
		
		k = unique(x$annotations[ 
      which(x$annotations[,'event'] %in% filter.in.names, arr.ind=T), 'event' 
      ])
		
		cat(paste(' [', length(k), '/', length(filter.in.names), ' found].', sep=''))

		valid[colnames] = TRUE
	}
	
	if(!any(is.na(filter.out.names)))
	{
		cat(paste('\nEvents for the following genes won\'t be selected (filter.out.names): ',
              sep='', paste(filter.out.names, collapse=', ')))

		colnames = which(x$annotations[,2] %in% filter.out.names, arr.ind=T)
		cat(paste(' [', length(colnames), '/', length(filter.out.names), ' found].', sep=''))
		
		valid[colnames] = FALSE
	}

	y = list()
	y$genotypes = x$genotypes[, valid]	
	
	y$annotations = as.matrix(x$annotations[valid, ])
	colnames(y$annotations) = c('type', 'event')
	rownames(y$annotations) = colnames(y$genotypes)

	y$types = as.matrix(x$types[unique(y$annotations[,1]), 1])
	colnames(y$types) = c('color')
	rownames(y$types) = unique(y$annotations[,1])

	if(!is.null(x$stages)) y$stages=x$stages
	is.compliant(x, err.fun='events.selection: output')
	
	cat(paste('\nSelected ', ncol(y$genotypes), ' events, returning.\n', sep=''))

	return(y)
}