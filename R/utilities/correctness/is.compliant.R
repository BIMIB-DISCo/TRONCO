# Check if 'x' is compliant with TRONCO's input: that is if it has dataframes x$genotypes, x$annotations, x$types and x$stage (optional)
# - x: list of dataframe to check
# - err.fun: string which identifies the function which called is.compliant
# - stage: boolean flag to check x$stage datagframe
#
# Returns: on error stops the computation
is.compliant = function(x, err.fun, stage=has.stages(x))
{
	# Check if x is defined
  	if(is.null(x) || is.na(x))
		stop(paste(err.fun, ': input \'x\' is null.'))

	# Check then if x is a list
  	if(!is.list(x))
		stop(paste(err.fun, ': input \'x\' is not a list.'))

	# Check then if x has required fields
    if(is.null(x$genotypes) || is.na(x$genotypes))
		stop(paste(err.fun, ': input \'x\' has no genotypes field.'))	
	else
		if(!is.matrix(x$genotypes) && !is.data.frame(x$genotypes))	stop(paste(err.fun, ': attribute genotypes in  \'x\' is not a matrix.'))	
	
	colnames(x$annotations) = c('type', 'event')
	colnames(x$types) = c('color')
	
	if(is.null(x$annotations) || is.na(x$annotations))
		stop(paste(err.fun, ': input \'x\' has no annotations field.'))
	else
		if(!is.matrix(x$annotations) && !is.data.frame(x$annotations)) stop(paste(err.fun, ': attribute annotations in  \'x\' is not a matrix.'))	

	if(is.null(x$types) || is.na(x$types))
		stop(paste(err.fun, ': input \'x\' has no types field.'))
	else
		if(!is.matrix(x$types) && !is.data.frame(x$types))  stop(paste(err.fun, ': attribute types in  \'x\' is not a matrix.'))	
	
	if(stage == TRUE && (is.null(x$stage) || is.na(x$stage)))
		stop(paste(err.fun, ': input \'x\' has no stage field.'))
	else
		if(stage == TRUE && !is.matrix(x$stage) && !is.data.frame(x$stage)) stop(paste(err.fun, ': attribute stage in  \'x\' is not a matrix.'))	

	# Annotations sould be present for all genotypes columns
  	if(nrow(x$annotations) != ncol(x$genotypes) ) 
		stop(paste(err.fun, ': input \'x\' has less annotations than expected.'))			
 	if(!all(colnames(x$genotypes) == rownames(x$annotations))) 
		stop(paste(err.fun, ': input \'x\' has inconsistent annotations.'))			

	# Types should be defined for every annotation
   	if(nrow(x$types) != length(unique(x$annotations[,1]))) 
		stop(paste(err.fun, ': input \'x\' has less types than expected.'))			

 	if(!all(unique(x$annotations[,'type']) %in% rownames(x$types))) 
	{
		stop(paste(err.fun, ': input \'x\' has inconsistent types (', 
				paste(unique(x$annotations[,'type']), collapse=',') 
				,' vs ', paste(rownames(x$types), collapse=',') ,').', 
			sep=''))				
	}
		
 	# Stage should be defined for every samples
  if(stage == TRUE && nrow(x$stages) != nrow(x$genotypes)) 
		stop(paste(err.fun, ': input \'x\' has less stages than expected.'))			
 	if(stage == TRUE && !all(rownames(x$stages) == rownames(x$genotypes))) 
		stop(paste(err.fun, ': input \'x\' has inconsistent stages.'))			
		
	if(stage == TRUE)	colnames(x$stages) = c('stage')
  
  if(has.duplicates(x)) {
    print("Duplicated events:")
    print(duplicates(x))
    warning('duplicated events found in annotations')
  }
 }