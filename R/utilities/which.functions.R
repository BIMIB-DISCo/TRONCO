which.samples = function(x, gene, type, neg = FALSE)
{
	data = as.gene(x, genes = gene, types = type)
	data = data[data == 1, , drop = FALSE] 
	
	samples = as.samples(x)
	
	if(neg)	return(setdiff(samples, rownames(data)))
	else return(rownames(data))	
}