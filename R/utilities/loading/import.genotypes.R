# Transforms a matrix 'geno' in an input list conformant to TRONCO's specifications. If this casting is not possible
# errors are thrown. If no geno.annot/stage.annot parameters are defined, column names of geno are used as event names.
# If these parameters are specified, they should be compliant with TRONCO's input (explain this better)
#
# - geno: a dataframe with the genotypes (this is not constrained to be 0/1)
# - stage.annot: a nrow(geno)x2 dataframe where stage.annot[i,] denotes the stage associated to geno[i,]. 
# - event.type: type label to be used
# - color: R Brewer palette from which random colors are sampled
# Returns: a list of dataframes compliant with TRONCO's specifications.
import.genotypes = function(geno, stage.annot = NA, event.type = 'variant', color = 'Accent')
{
	if (!require('RColorBrewer')) {
    	install.packages('RColorBrewer', dependencies = TRUE)
    	library(RColorBrewer)
  	}
  	
	# Avoid malformed datasets
  	if(ncol(geno) == 0 || nrow(geno) == 0) stop('Empty genotypes (number of rows/columns 0), will not import.')
	nc = ncol(geno)
	nr = nrow(geno)
	
	# Gather col/row names
	if(is.null(colnames(geno))) 
	{
		cn = paste0('Gene', 1:ncol(geno))
		warning('Missing column names to identify genes. Will use labels \"Gene1\", \"Gene2\", .....')
	} else 	cn = colnames(geno)

	# Gather col/row names
	if(is.null(rownames(geno))) 
	{
		rn = paste0('Sample', 1:nrow(geno))
		warning('Missing row names to identify samples. Will use labels \"Sample1\", \"Sample2\",  .....')
	} else 	rn = rownames(geno)

	x = list()
	
	# Access keys - G1, G2, ...
	keys = paste0('G', 1:ncol(geno))
	
	# Genotype matrix
	x$genotypes = as.matrix(geno)
	colnames(x$genotypes) = keys
	rownames(x$genotypes) = rn
	
	# Create attributes
	x$annotations = matrix(0, nrow=nc, ncol=2)
	colnames(x$annotations) = c('type', 'event')
	rownames(x$annotations) = keys
	
	x$annotations[, 'type']  = event.type
	x$annotations[, 'event'] = cn
	
	# We create a map from types to colors
	x$types = matrix(color, nrow=1, ncol=1)
	rownames(x$types) = event.type
	colnames(x$types) = c('color')
	
	# If the input color is a ColorBrewer scheme	
	#	my.palette = color	
	#if(color %in% rownames(brewer.pal.info)) 
	#	my.palette = brewer.pal(n=brewer.pal.info[color, 'maxcolors'], name=color); 
	#
			
	is.compliant(x, 'import.genotypes: output')		

	return(x)
 }