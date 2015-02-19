#
# import.gistic - Convert a GISTIC score file to TRONCO input.
#

# TODO - Stages?

"import.gistic" <- function(x, merge = FALSE, amp.del = FALSE) {

	cat('*** GISTIC input format conversion started.\n'); 	

	if(merge == TRUE && merge == amp.del) 
		stop('merge = amp.del = TRUE')

	# For next operations it is convenient to have everything as 'char' rather than 'int'
	# if(typeof(x) == typeof(list()))
	# {
		# rn = rownames(x)
		# x= apply(x, 2, as.character) 
		# rownames(x) = rn		
	# }	
	 # if(typeof(x) == typeof(list()))
		# cat('Warning: maybe convert input to char should speed-up computation.')
	
	# x <- data.frame(lapply(x, as.integer), stringsAsFactors=FALSE, row.names=rownames(x))
    # x = as.matrix(sapply(x, as.numeric), stringsAsFactors=FALSE, row.names=rownames(x))  

	if(any(is.na(x))) cat('Some entries are either NA or NaN and will be replaced with 0s.\n')
	x[is.na(x)] = 0		
    
	k <- 2
	if(!merge) k = 4;
    cat(paste('#events to define (at most): ', k * (ncol(x)), '\n', sep=''));
	
	# gene symbols
	enames <- colnames(x)
	if(is.null(enames)) stop('Error: gistic file has no column names and can not imported, aborting!')
	
	cat(paste('Gene Symbols: ', sep='')); 
	cat(str(enames))

	cat('\nExtracting: Homozygous Loss, '); 	
	d.homo = x
	d.homo[d.homo != -2] <- 0
	d.homo[d.homo == -2] <- 1
		
	cat('Heterozygous Loss, '); 	
	d.het  <- x
	d.het[d.het != -1] <- 0
	d.het[d.het == -1] <- 1
		
	cat('Low-level Gain and '); 	
	d.low  <- x
	d.low[d.low != 1] <- 0

	cat('High-Level Gain.\n'); 	
	d.high <- x
	d.high[d.high != 2] <- 0
	d.high[d.high == 2] <- 1
	
	if(amp.del)
	{
		cat('Amplifications (High-Level Gain) and Deletions (Homozygous Loss) will be returned.\n')
		
		ampl = import.genotypes(d.high, default.variant='Amplification', color = 'firebrick4')
		del = import.genotypes(d.homo, default.variant='Deletion', color = 'dodgerblue4')
		
		d.cnv.all = ebind(ampl, del)
		
		is.compliant(d.cnv.all, 'import.gistic: output')
		return(d.cnv.all);
	}
	
	if(merge)	
	{
		# This is the 'slow' solution 
		# .... = union.types(d.cnv.all, 'Homozygous Loss', 'Heterozygous Loss', new.type='Loss')
		# .... = union.types(d.cnv.all, 'Low-level Gain', 'High-level Gain', new.type='Gain')	
	
		cat('Homozygous and Heterozygous loss will be merged. ')
		cat('Low-level and High-Level gain will be merged.\n')
		
		# This is quicker: convert everething to int (could be char)
		rn = rownames(d.homo)
		cn = colnames(d.homo)
		d.homo = apply(d.homo, 2, as.integer) 
		d.het = apply(d.het, 2, as.integer) 
		d.low = apply(d.low, 2, as.integer) 
		d.high = apply(d.high, 2, as.integer) 
		rownames(d.homo) = rn
		colnames(d.homo) = cn		
		 	
		# Implement OR via sum 	
		loss = d.homo + d.het
		loss[loss == 2] = 1
		
		gain = d.low + d.high
		gain[gain == 2] = 1
		
		loss = import.genotypes(loss, default.variant='Loss', color='dodgerblue4')
		gain = import.genotypes(gain, default.variant='Gain', color='firebrick4')
		
		d.cnv.all = ebind(loss, gain)
	}
	else
	{
		homo = import.genotypes(d.homo, default.variant='Homozygous Loss')
		homo$types['Homozygous Loss', ] = 'dodgerblue4'
		
		het = import.genotypes(d.het, default.variant='Heterozygous Loss')
		het$types['Heterozygous Loss', ] = 'dodgerblue1'
		
		low = import.genotypes(d.low, default.variant='Low-level Gain')
		low$types['Low-level Gain', ] = 'firebrick1'

		high = import.genotypes(d.high, default.variant='High-level Gain')
		high$types['High-level Gain', ] = 'firebrick4'
	
		d.cnv.all = ebind(homo, het, low, high)
	}
	
	cat(paste('*** Data extracted, returning ', ncol(d.cnv.all$genotypes),
		' events (', length(unique(d.cnv.all$annotations[, 'event'])),' genes) and ', 
		nrow(d.cnv.all$genotypes), ' samples.\n', sep=''));
		
	is.compliant(d.cnv.all, 'import.gistic: output')
	# d.cnv.all <- d.cnv.all[,which(!apply(d.cnv.all ==0,2,all))]
	# print('These events were found in at least one patient')
	# print(colnames(d.cnv.all))	
    # print(paste('Boolean matrix size: ', nrow(d.cnv.all), ' x ', ncol(d.cnv.all), sep=''));

	return(d.cnv.all);
}
