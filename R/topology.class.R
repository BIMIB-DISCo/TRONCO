##################################################################################
#                                                                                #
# TRONCO: a tool for TRanslational ONCOlogy                                      #
#                                                                                #
##################################################################################
# Copyright (c) 2014, Marco Antoniotti, Giulio Caravagna, Alex Graudenzi,        #
# Ilya Korsunsky, Mattia Longoni, Loes Olde Loohuis, Giancarlo Mauri, Bud Mishra #
# and Daniele Ramazzotti.                                                        #
#                                                                                #
# All rights reserved. This program and the accompanying materials               #
# are made available under the terms of the Eclipse Public License v1.0          #
# which accompanies this distribution, and is available at                       #
# http://www.eclipse.org/legal/epl-v10.html and in the include COPYING file      #
#                                                                                #
# Initial contributors:                                                          #
# Giulio Caravagna, Alex Graudenzi, Mattia Longoni and Daniele Ramazzotti.       #
##################################################################################

#' @import methods
# The topology class definition
topology = NULL;
setClass("topology",
		representation(			  
						dataset = "data.frame",
                        
					   	marginal.probs = "matrix",
					   	joint.probs = "matrix",
					    cond.probs = "matrix",
                       
				 		estimated.marginal.probs = "matrix",
					    estimated.joint.probs = "matrix",
					    estimated.cond.probs = "matrix",
                        edge.confidence = "matrix",
                        confidence = "list",
         
					    pr.score = "matrix",
					    adj.matrix = "matrix",
					    adj.matrix.bic = "matrix",
            
					    is.valid = "logical",
                        invalid.events = "list",
                       
					    error.fp = "numeric",
                        error.fn = "numeric",
                        bootstrap.settings = "list",
                        bootstrap = "logical",
                        
                        algorithm = "character"))

# A summary is displayed if the object name of class 
# topology is written in the R console.
setMethod("show", "topology",
		  function(object){
		  	tree <- c()
		  	adj.matrix <- object@adj.matrix
		  	names <- colnames(adj.matrix)
		  	for(i in 1:nrow(adj.matrix))
		  		for(j in 1:ncol(adj.matrix))
		  			if(adj.matrix[i,j] == 1)
		  				tree <- c(tree, 
		  						  paste(names[i], " -> ", names[j], "\n"))
		  	cat(" Tree progression model with ", ncol(adj.matrix), "events \n")
		  	cat(" ")
		  	cat(tree)
		  	cat("\n Estimated false positives error rate:", object@error.fp)
		  	cat("\n Estimated false negative error rate:", object@error.fn)
		  	
        if(object@bootstrap){
          cat("\n")
          cat("\n Executed", object@bootstrap.settings$type, "bootstrap, nboot:", topology@bootstrap.settings$nboot)
          cat("\n Confidence overall value:", object@confidence$overall.value)
          cat("\n Confidence overall frequency:", object@confidence$overall.frequency)
        }
})