##################################################################################
#                                                                                #
# TRONCO: a tool for TRanslational ONCOlogy                                      #
#                                                                                #
##################################################################################
# Copyright (c) 2015, Marco Antoniotti, Giulio Caravagna, Luca De Sano,          #
# Alex Graudenzi, Ilya Korsunsky, Mattia Longoni, Loes Olde Loohuis,             #
# Giancarlo Mauri, Bud Mishra and Daniele Ramazzotti.                            #
#                                                                                #
# All rights reserved. This program and the accompanying materials               #
# are made available under the terms of the GNU GPL v3.0                         #
# which accompanies this distribution                                            #
#                                                                                #
##################################################################################


# Check if 'x' is compliant with TRONCO's input: that is if it has dataframes x$genotypes, x$annotations, x$types and x$stage (optional)
# - x: list of dataframe to check
# - err.fun: string which identifies the function which called is.compliant
# - stage: boolean flag to check x$stage datagframe
#
# Returns: on error stops the computation
#' @export
is.compliant = function(x, err.fun='[ERR]', stage = !(all(is.null(x$stages)) || all(is.na(x$stages))))
{
	# Check if x is defined
  	if(is.null(x) || is.na(x))
		stop(paste(err.fun, ': input \'x\' is null.'))

	# Check then if x is a list
  	if(!is.list(x))
		stop(paste(err.fun, ': input \'x\' is not a list.'))

	# Check then if x has required fields
    if(is.null(x$genotypes) || any(is.na(x$genotypes)))
		stop(paste(err.fun, ': input \'x\' has no genotypes field.'))	
	else
		if(!is.matrix(x$genotypes) && !is.data.frame(x$genotypes))	stop(paste(err.fun, ': attribute genotypes in  \'x\' is not a matrix.'))	
	
	colnames(x$annotations) = c('type', 'event')
	colnames(x$types) = c('color')
	
	if(is.null(x$annotations) || any(is.na(x$annotations)))
		stop(paste(err.fun, ': input \'x\' has no annotations field.'))
	else
		if(!is.matrix(x$annotations) && !is.data.frame(x$annotations)) stop(paste(err.fun, ': attribute annotations in  \'x\' is not a matrix.'))	

	if(is.null(x$types) || any(is.na(x$types)))
		stop(paste(err.fun, ': input \'x\' has no types field.'))
	else
		if(!is.matrix(x$types) && !is.data.frame(x$types))  stop(paste(err.fun, ': attribute types in  \'x\' is not a matrix.'))	
	
	if(stage == TRUE && (is.null(x$stages) || all(is.na(x$stages))))
		stop(paste(err.fun, ': input \'x\' has no stage field.'))
	else
		if(stage == TRUE && !is.matrix(x$stages) && !is.data.frame(x$stages)) stop(paste(err.fun, ': attribute stage in  \'x\' is not a matrix.'))	

	# Annotations sould be present for all genotypes columns
  	if(nrow(x$annotations) != ncol(x$genotypes) ) 
		stop(paste(err.fun, ': input \'x\' has less annotations than expected.'))			
 	if(!all(colnames(x$genotypes) == rownames(x$annotations))) 
		stop(paste(err.fun, ': input \'x\' has inconsistent annotations.'))			

	# Types should be defined for every annotation
   	if(nrow(x$types) != length(unique(x$annotations[,1]))) 
	{
		message('[ERROR] rownames(x$types):', paste(rownames(x$types), collapse=', ', sep=''))
		message('[ERROR] Annotated event types:', paste(unique(x$annotations[,1]), collapse=', ', sep=''))
	
		stop(paste(err.fun, ': input \'x\' has less types than expected.'))			
	}
	
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
  
  dup = duplicated(x$annotations)
  if(any(dup)) 
    {
    cat("Duplicated events in \'x\': \n")
    print(head(x$annotations[dup, ]))
    stop('Duplicated events.')
  }
 }


# Check if x is a valid TRONCO model
is.model = function(x)
{
  if(!'model' %in% names(x))
    stop('Input object is not a TRONCO model.')
}


# Check if y is a valid event list for x
is.events.list = function(x, y)
{
  if(!is.matrix(y)) stop('Events should be given as matrix - see "as.events".')
  if(ncol(y) != 2 ||
       !all(c('type', 'event') %in% colnames(y))
       ) stop('Events are missing column "type" (type of event) or "event" (gene symbol) - see "as.events".')
  
   if (!all(rownames(y) %in% colnames(x$genotypes))) 
    stop('Events rownames are not valid keys for genotypes - see "as.events".')
}
