#
# import.mutations - Format supported is ....
#

"import.mutations" <- function(x, merge = FALSE, color=NA) {

	cat('*** Mutations input format conversion started.\n'); 	
	x <- data.frame(lapply(x, as.character), stringsAsFactors=FALSE, row.names=rownames(x))

	if(any(is.na(x))) cat('Found NA/NaN entries which and will be replaced with 0s.\n')
	x[is.na(x)] = 0		
    x[x==NaN] = 0
    
    if(!merge) x[ x != 0 ] = 1
    
    if(is.na(color)) x = import.genotypes(data.matrix(x), default.variant='SNV')
    else x = import.genotypes(data.matrix(x), default.variant='SNV', color=color)
    
	cat(paste('*** Data extracted, returning ', ncol(x$genotypes),' events and ', nrow(x$genotypes), ' samples.\n', sep=''));
	is.compliant(x)

	return(x);
}
