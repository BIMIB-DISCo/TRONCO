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

library(lattice)
#' @import lattice
#' @export tronco.bootstrap.show
#' @title show bootstrapping results
#'
#' @description
#' \code{tronco.bootstrap.show} show bootstrapping results. Requires that you already executed tronco.bootstrap 
#' 
#' @param topology A topology returned by a reconstruction algorithm
tronco.bootstrap.show <- function(topology){
  if(missing(topology))
    stop("Missing parameter for tronco.bootstrap.show function: tronco.bootstrap.show(topology)", call. = FALSE)
	if(!topology@bootstrap)
      stop("To show confidence information bootstrap execution is needed! see: tronco.bootstrap function!", call. = FALSE)
	print(topology@edge.confidence)
	levelplot(topology@edge.confidence, xlab = "", ylab = "", 
	scales = list(x = list(alternating = 2, rot = 90), tck = 0), 
	main = paste("Edge confidence (", topology@bootstrap.settings$type, " bootstrap)",sep = ""))
}