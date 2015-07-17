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

#' @export
which.samples = function(x, gene, type, neg = FALSE)
{
	data = as.gene(x, genes = gene, types = type)
	data = data[data == 1, , drop = FALSE] 
	
	samples = as.samples(x)
	
	if(neg)	return(setdiff(samples, rownames(data)))
	else return(rownames(data))	
}
