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
	cat(paste('\nFilters : freq. [', !is.na(filter.freq), '] \t names.in [', !any(is.na(filter.in.names)), ']\t names.out [',
            !any(is.na(filter.out.names)), ']', sep=''))

	if(is.na(filter.out.names) && is.na(filter.in.names) && is.na(filter.freq))
		return(x);

	valid = rep(FALSE, ncol(x$genotypes))
					
	if(!is.na(filter.freq))
	{	
		cat(paste('\nExtracting events with a minimum frequency of ', filter.freq, ' (', 
              round(nsamples(x) * filter.freq, 0),' alterations out of ', nsamples(x),' samples).\n', sep=''))
    x = enforce.numeric(x)		
    
		pb = txtProgressBar(1, nevents(x), style = 3);      
		for(i in 1:nevents(x))
		{		
		  setTxtProgressBar(pb, i)  
		  
			mut.freq <- sum(x$genotypes[,i])/nsamples(x)
			valid[i] <- mut.freq > filter.freq
		}
		close(pb)

    nev = min(10, nrow(as.events(x)[valid, ]))
    cat(paste('Selected ', nrow(as.events(x)[valid, ]), ' events, showing ', nev, '.\n', sep=''))
    print(head(as.events(x)[valid, ], n=nev))
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
	
	cat(paste('\nSelected ', nevents(y), ' events, returning.\n', sep=''))

	return(y)
}

# Return the first n recurrent events
rank.recurrents = function(x, n)
{
	is.compliant(x)
    x = enforce.numeric(x)		

	
	if(n <= 0) stop('Rank value (n) should be positive.')
	
	# Sum columns
	sums = colSums(x$genotypes)
	
	# Get the names of the first n ranked
	sorted = sort(sums, decreasing = T)
	
	# print(sorted[1:20])
	
	scores = unique(sorted)
	# print(scores)	

	l = length(scores)
	if(n >l) warning(paste0('Rank contains ', l, ' unique entries, using n=', l, ' instead of n=', n))

	n = min(n, length(scores))
	scores = scores[1:n]
	
	sorted = sorted[which(sorted >= min(scores))]

	max = names(sorted[which(sorted == max(scores))])
	min = names(sorted[which(sorted == min(scores))])

	cat(paste0('Most recurrent(s): ', paste(as.events(x)[max, 'event'], collapse=', '), ' (', (max(scores)), ' hits).\n' ))
	cat(paste0(n, '-th recurrent(s): ', paste(as.events(x)[min, 'event'], collapse=', '), ' (', (min(scores)), ' hits).\n' ))
	
	
	order = names(sorted)
	genes = as.events(x)[order, 'event']
	
	return(as.vector(genes))
	
}