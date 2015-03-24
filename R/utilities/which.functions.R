which.samples = function(x, gene, type, neg = FALSE)
{
	data = as.genotypes(x)
	
	keys = rownames(as.events(x, genes=event, types = type))
	
	data = data[, keys, drop = FALSE]
	samples = as.samples(x)
	data = data[data == 1, ] 
	
	if(neg)	return(setdiff(samples, rownames(data)))
	else return(rownames(data))	
}