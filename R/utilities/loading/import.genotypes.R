# Transforms a matrix 'geno' in an input list conformant to TRONCO's specifications. If this casting is not possible
# errors are thrown. If no geno.annot/stage.annot parameters are defined, column names of geno are used as event names.
# If these parameters are specified, they should be compliant with TRONCO's input (explain this better)
#
# - geno: a dataframe with the genotypes (this is not constrained to be 0/1)
# - geno.annot: a ncol(geno)x2 dataframe where geno.annot[i,] denotes the event associated to geno[,i]. Column should contain event type and label
# - stage.annot: a nrow(geno)x2 dataframe where stage.annot[i,] denotes the stage associated to geno[i,]. 
# - default.variant: default type label used if geno.annot=NA
# - color: R Brewer palette from which random colors are sampled
# Returns: a list of dataframes compliant with TRONCO's specifications.
import.genotypes = function(geno, geno.annot=NA, stage.annot=NA, default.variant='var.', color='Accent')
{
	if (!require('RColorBrewer')) {
    	install.packages('RColorBrewer', dependencies = TRUE)
    	library(RColorBrewer)
  	}

	# Check if the input is a dataframe or a file to read from disk
	
	
	if(ncol(geno) == 0 || nrow(geno) == 0)
		stop('Genotypes have wrong number of rows/columns (0).')
	
	x = list()
	
	# Copy genotype matrix, and sets its rownames (samples)
	x$genotypes = as.matrix(geno)
	rownames(x$genotypes) = rownames(geno)	
		
	# Names are stored as attributes
	xn = colnames(geno)
	nc = ncol(geno)
	nr = nrow(geno)
	
	# Create attributes
	x$annotations = matrix(0, nrow=nc, ncol=2)
	colnames(x$annotations) = c('type', 'event')
	
	if(is.na(geno.annot))
	{
		# If no attribute is given, default 'variant' is used
		x$annotations[, 'type']  = default.variant
		x$annotations[, 'event'] = xn
	}
	else # Otherwise, use input
		x$annotations = geno.annot # TODO: debug this...
					
	# We create a mapping from columns in x$genotypes and row names in x$attributes
	names.map = paste('G', c(1:nc))
	colnames(x$genotypes) = names.map
	rownames(x$annotations) = names.map
	
	# We create a map from types to colors
	num.types = length(unique(x$annotations[, 'type']))
	
	x$types = matrix(0, nrow=num.types, ncol=1)
	rownames(x$types) = unique(x$annotations[,1])
	colnames(x$types) = c('color')

	# If the input color is a ColorBrewer scheme
	my.palette = color;
	
	if(color %in% rownames(brewer.pal.info)) 
		my.palette = brewer.pal(n=brewer.pal.info[color, 'maxcolors'], name=color); 

	x$types[, 'color'] = colorRampPalette(my.palette)(num.types)
	
	
	# Add stage, if given as input
	if(!any(is.na(stage.annot)))
	{
		if(nrow(stage.annot) != nr)
			cat(paste('Missing stage information for some samples (genotypes= ', nr, ', stages=', nrow(stage.annot),'), setting them as NA.\n', sep=''))
				
		x$stages = matrix('NA', nrow=nr, ncol=1)	
		# x$stage = data.frame(row.names=rownames(x$genotypes))
		# x$stage = cbind(as.character(stage.annot[rownames(x$stage), 1]))		
		# print(x$stage)
		# print(stage.annot[rownames(x$genotypes), 1])
		rownames(x$stages) = rownames(x$genotypes)
		# print(x$stage)
		
		# print(stage.annot)
		x$stages[,1] = as.character(stage.annot[rownames(x$genotypes), 1])
		# print(x$stage)
		# print(rownames(x$stage))
		
		rownames(x$stages) = rownames(x$genotypes)
	}
	# print(x)
	# else x$stage = NULL
	# print(nrow(x$genotypes))
	# print(nrow(as.matrix(x$stages)))
	
	# print('sss')
	# print(as.matrix(x$stage))
	 # print(nrow(x$stages))
	 # print(nrow(x$genotypes))
	 # print(dim(x$genotypes))
	
	 # print(nrow(x$stages) != nrow(x$genotypes))
	
	# print(is.na(stage.annot))
	# print(is.null(stage.annot))
	
	if(is.na(stage.annot))
		is.compliant(x, 'import.genotypes', stage=FALSE)		
	else
    	is.compliant(x, 'import.genotypes', stage=TRUE)		
	
	return(x)
 }