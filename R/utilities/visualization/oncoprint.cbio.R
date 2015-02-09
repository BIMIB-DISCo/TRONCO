#
# oncoprint.cbio : export input for cbio visualization at http://www.cbioportal.org/public-portal/oncoprinter.jsp
#
oncoprint.cbio <- function(x, file='oncoprint-cbio.txt', hom.del = 'Homozygous Loss',het.loss = 'Heterozygous Loss', gain = 'Low-level Gain', amp = 'High-level Gain')
{
	is.compliant(x)
	
	# r = paste(paste(rownames(x$genotypes)), x$annotations[x$annotations[,''],], 'xxx')
	r = 'Sample\tGene\tAlteration\n'
	for(i in 1:nrow(x$genotypes))
	{
		for(j in 1:ncol(x$genotypes))
		{
			if(x$genotypes[i,j] == 1)
			{
				s = rownames(x$genotypes)[i]
				g = x$annotations[colnames(x$genotypes)[j], 'event']
				
				t = x$annotations[colnames(x$genotypes)[j], 'type']
				
				t.o = 'xxx'
				if(t == hom.del) t.o =  'HOMDEL'
				if(t == het.loss)  t.o =  'HETLOSS'
				if(t == gain)  t.o =  'GAIN'
				if(t == amp)  t.o =  'AMP'
				
				# cat(paste( s,  g,  t.o, '\n', sep=' ', collpase=''))
				r = paste(r, s, '\t', g, '\t', t.o, '\n', sep='', collpase='')
			}
		}	
		
	}
	
	write(r, file)
}

