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

#' @export tronco.caprese
#' @title runs CAPRESE algorithm
#'
#' @description
#' \code{tronco.caprese} executes the CAPRESE algorithm on the dataset \code{data.values} specified. 
#'
#' @details
#' \code{tronco.caprese} executes the reconstruction of the topology, and computesg all the confidence measures defined in \code{confidence}.
#' 
#' @param dataset The input dataset. Type: dataframe. The dataset given as input is the data.values data frame loaded by the \code{data} function.
#' @seealso \code{\link{data}}
#' @param lambda the real positive value of the shrinkage coefficient, required to range in [0, 1]. Its default value is 0.5, if unspecified.
#' @param verbose execute CAPRESE algorithm with verbose output to screen. Type: boolean, dafault: FALSE.
#' @return an object containing the reconstructed topology and confidence values.
#' 
tronco.caprese <- function(dataset, lambda = 0.5, verbose = FALSE){
  
  if(missing(dataset))
    stop("Missing parameter for tronco.caprese function: tronco.caprese(dataset, lambda, verbose)", call. = FALSE)
  if(lambda < 0 || lambda > 1)
    stop("Lambda coefficient must be in [0:1]!", call. = FALSE)
  else{
    if(is.null(dataset))
      stop("Empty dataset!", call. = FALSE)
    
    # Load CAPRESE algorithm
    sapply(list.files(pattern="[.]R$", path="CAPRESE", full.names=TRUE), source)
    
    # Reconstruct the topology
    topology <- caprese.fit(dataset, lambda, verbose)

    
    	if(!exists("invalid.events"))
    		stop("Invalid.events collection not fount! use data.load() function to load a dataset and looks for invalid events...", call = FALSE)
       invalid.events <- invalid.events
       # Collects info for the object of class topology.
       info  <- infoset(invalid.events$merged.events, invalid.events$removed.events)
       
       rownames(topology$adj.matrix) <- info$all.labels
       colnames(topology$adj.matrix) <- info$all.labels
       
       rownames(topology$pr.score) <- info$all.labels
       colnames(topology$pr.score) <- info$all.labels
       
       # Probability Labels assignment
       rownames(topology$probabilities$conditional.probs) <- info$all.labels
       rownames(topology$probabilities$marginal.probs) <- info$all.labels
       
       rownames(topology$probabilities$joint.probs) <- info$all.labels
       colnames(topology$probabilities$joint.probs) <- info$all.labels
       
       # Estimated probability labels assignment
       rownames(topology$probabilities$estimated.conditional.probs) <- info$all.labels
       rownames(topology$probabilities$estimated.marginal.probs) <- info$all.labels
       
       rownames(topology$probabilities$estimated.joint.probs) <- info$all.labels
       colnames(topology$probabilities$estimated.joint.probs) <- info$all.labels
       
       
       topology.obj <- new("topology",
    	   	 					   dataset = topology$dataset,
                         
    	   	 					   marginal.probs = topology$probabilities$marginal.probs,
    	   	 					   joint.probs = topology$probabilities$joint.probs,
    	   	 					   cond.probs = topology$probabilities$conditional.probs,
                       
    	   	 					   estimated.marginal.probs = topology$probabilities$estimated.marginal.probs,
    	   	 					   estimated.joint.probs = topology$probabilities$estimated.joint.probs,
    	   	 					   estimated.cond.probs = topology$probabilities$estimated.conditional.probs,
                       
    	   	 					   edge.confidence = matrix(),
    	   	 					   confidence = list(),
    	   	 					   bootstrap.settings = list(),
    	   	 					   bootstrap = FALSE,
                         
    	   	 					   pr.score = topology$pr.score,
    	   	 					   adj.matrix = topology$adj.matrix,
                       		   adj.matrix.bic = matrix(),
                       		   
                         			# The parameter is fixed because of the last code changes ;)
    	   	 					   is.valid = TRUE,
                                invalid.events = invalid.events,
                         
    	   	 					   error.fp = topology$error.rates$error.fp,
    	   	 					   error.fn = topology$error.rates$error.fn,
                       		   algorithm = "CAPRESE")
       
       cat(paste("Executed CAPRESE algorithm with shrinkage coefficient:", lambda,"\n"))
       cat(" Estimated false positives error rate: ", topology.obj@error.fp)
       cat("\n Estimated false negative error rate: ", topology.obj@error.fn)
       
       return(topology.obj)
       
    
  }
  
}