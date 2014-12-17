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

#' @import lattice
#' @export tronco.bootstrap
#' @title perform bootstrap algorithm
#'
#' @description
#' \code{tronco.bootstrap} perform parametric and non-parametric bootstrap algorithms 
#' 
#' @param topology A topology returned by a reconstruction algorithm
#' @param lambda A lambda value, default is 0.5
#' @param type The type of bootstrap performed, parametric and non parametric types are available.
#' To specify wich type of bootstrap run type must be "parametric" or "non-parametric".
#' @param nboot Samplig value. The grater will be the nboot value the logehr time the
#' entire process will take to complete the computing
#' @return
#' A topology object with bootstrap informations added
#' 
tronco.bootstrap <- function(topology, lambda = 0.5, type = "non-parametric", nboot = 1000){
   show.level <- function(topology){
		print(topology@edge.confidence)
		levelplot(topology@edge.confidence, xlab = "", ylab = "", 
			scales = list(x = list(alternating = 2, rot = 90), tck = 0), 
			main = paste("Edge confidence (", topology@bootstrap.settings$type, " bootstrap)",sep = ""))
   }
   
   if(missing(topology))
     stop("Missing parameter for tronco.bootstrap function: tronco.bootstrap(topology, lambda, type, nboot)", call. = FALSE)
   # If the given topology is valid and contains a valid dataset and other informations.
   if(topology@is.valid){
   	 # If all numeric parameters are well formed.
     if(nboot > 0){
       if(lambda > 0 && lambda < 1){
         if(type == "non-parametric" || type == "parametric"){
           	
           cat("Executing bootstrap algorithm this may take several time...\n")
           
           # The error rates are turned into the format required by bootstrap.caprese function.
           error.rates <- list()
           error.rates$error.fp <- topology@error.fp
           error.rates$error.fn <- topology@error.fn
           
           # The bootstrap function of the caprese algorithm is performed.
           boot.info <- bootstrap.caprese(topology@dataset, lambda, topology@adj.matrix, type, topology@estimated.marginal.probs, topology@estimated.cond.probs, error.rates, nboot)
           
           # Confidence informations are set into the topology object.
           topology@edge.confidence <- boot.info$edge.confidence
           topology@confidence <- boot.info$confidence
           topology@bootstrap.settings <- boot.info$bootstrap.settings
           topology@bootstrap <-  TRUE
           
           cat(paste("Executed ", type, " bootsrap with ", nboot, 
                     " as sampling number and ", lambda, " as lambda value\n\n", sep =""))
           
           print(topology@edge.confidence)
           cat("\nConfidence overall value:", boot.info$confidence$overall.value)
           cat("\nConfidence overall frequency:", boot.info$confidence$overall.frequency)
     		 
         }
         else
           stop("Valid type of bootstrap are: non-parametric or parametric, check your typing!", call. = FALSE)
       }
       else
         stop("Lambda must be in [0:1]!", call. = FALSE)
     }
     else
       stop("Bootstrap number bust be greater than zero!", call. = FALSE)
   }
   else
     stop("Topology object does not contain a valid CAPRESE reconstruction!", call. = FALSE)
      
   return(topology)
}