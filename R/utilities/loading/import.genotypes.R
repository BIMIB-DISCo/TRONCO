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
	
	
					
	# We create a mapping from columns in x$genotypes and row names in x$attributes
	names.map = paste0('G', c(1:nc), sep='')
	colnames(x$genotypes) = names.map
	rownames(x$annotations) = names.map
	
	# We create a map from types to colors
	num.types = length(unique(x$annotations[, 'type']))
	
	x$types = matrix(0, nrow=num.types, ncol=1)
	rownames(x$types) = unique(x$annotations[, 1])
	colnames(x$types) = c('color')

	# If the input color is a ColorBrewer scheme
	my.palette = color;
	
	if(color %in% rownames(brewer.pal.info)) 
		my.palette = brewer.pal(n=brewer.pal.info[color, 'maxcolors'], name=color); 

	x$types[, 'color'] = colorRampPalette(my.palette)(num.types)
	
	
	# Add stage, if given as input
	if(!all(is.na(stage.annot)))
	{
		if(nrow(stage.annot) != nr)
			cat(paste('Missing stage information for some samples (genotypes= ', nr, ', stages=', nrow(stage.annot),'), setting them as NA.\n', sep=''))
				
		x$stages = matrix('NA', nrow=nr, ncol=1)	
		rownames(x$stages) = rownames(x$genotypes)
		x$stages[,1] = as.character(stage.annot[rownames(x$genotypes), 1])
		rownames(x$stages) = rownames(x$genotypes)
	}
		
	is.compliant(x, 'import.genotypes: output')		

	return(x)
 }