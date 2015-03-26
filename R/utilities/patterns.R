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
