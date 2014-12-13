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
#' @export confidence.data.joint
#' @export confidence.fit.joint
#' @export confidence.data.single
#' @export confidence.fit.single
#' @export confidence.data.conditional
#' @export confidence.fit.conditional
#' @export confidence.single
#' @export confidence.joint
#' @export confidence.conditional
#' @name confidence
#' @title provides various kinds of confidence measures for an inferred progression model
#' @description
#' A set of functions to visualise and compare the probability of each event in the progression model, as well as their joint and conditional distributions. These can be evaluated both in the data (observed probabilities) and in the reconstructed model (fitted probabilities).
#'
#' @param topology A topology returned by the reconstruction algorithm
#' @usage confidence.data.joint(topology) 
#' @details \code{confidence.data.joint} plot the pairwise observed joint probability of the events 
confidence.data.joint <- function(topology){
  
  if(missing(topology))
    stop("Missing parameter for confidence.data.joint function: confidence.data.joint(topology)", call. = FALSE)
	# Displays the joint.probs slot of the topology variable
	print(topology@joint.probs)
	# Plots a "levelplot" for joint.probs slot in topology variable
	levelplot(topology@joint.probs, xlab = "", ylab = "", scales = list(x = list(alternating = 2, rot = 90), tck = 0), main = "Data (joint probability)")
}

#' @rdname confidence
#' @usage confidence.fit.joint(topology) 
#' @details \code{confidence.fit.joint} plot the pairwise fitted joint probability of the events
confidence.fit.joint <- function(topology){
  if(missing(topology))
    stop("Missing parameter for confidence.fit.joint function: confidence.fit.joint(topology)", call. = FALSE)
  print(topology@estimated.joint.probs)
  levelplot(topology@joint.probs, xlab = "", ylab = "", scales = list(x = list(alternating = 2, rot = 90), tck = 0), main = "Fit (joint probability)")
}

#' @rdname confidence 
#' @usage confidence.data.single(topology)
#' @details \code{confidence.data.single} plot the observed probability of each event
confidence.data.single <- function(topology){
  if(missing(topology))
    stop("Missing parameter for confidence.data.single function: confidence.data.single(topology)", call. = FALSE)
	# Gets names from the marginal.prob matrix header
	names <- rownames(topology@marginal.probs) 
	v <- as.vector(topology@marginal.probs)
	# Sets the name for each element in vector v
	names(v) <- names
	# Sets borders for the plot window 
	par("mar" = c((max(nchar(names))/1.6), 2, 2, 1), "mfrow" = c(1,1))
	barplot(v, las = 3, main = "Data (events probability)")
	print(v)
}

#' @rdname confidence
#' @usage confidence.fit.single(topology)
#' @details \code{confidence.fit.single} plot the fitted probability of each event
confidence.fit.single <- function(topology){
  if(missing(topology))
    stop("Missing parameter for confidence.fit.single function: confidence.fit.single(topology)", call. = FALSE)
  names <- rownames(topology@estimated.marginal.probs)
  v <- as.vector(topology@estimated.marginal.probs)
  names(v) <- names
  #Elimino i bordi attorno al plot
  par("mar" = c((max(nchar(names))/1.6), 2, 2, 1), "mfrow" = c(1,1))
  barplot(v, las = 3, main = "Fit (events probability)")
  print(v)
}

#' @rdname confidence 
#' @usage confidence.data.conditional(topology)
#' @details \code{confidence.data.conditional} plot the pairwise observed conditional probability of the events
confidence.data.conditional <- function(topology){
  if(missing(topology))
    stop("Missing parameter for confidence.data.conditional function: confidence.data.conditional(topology)", call. = FALSE)
	names <- rownames(topology@cond.probs)
	v <- as.vector(topology@cond.probs)
	names(v) <- names
	par("mar" = c((max(nchar(names))/1.6), 2, 2, 1), "mfrow" = c(1,1))
	barplot(v, las = 3, main = "Data (conditional probability)")
	print(v)
}

#' @rdname confidence
#' @usage confidence.fit.conditional(topology)
#' @details \code{confidence.fit.conditional} plot the pairwise fitted conditional probability of the events
confidence.fit.conditional <- function(topology){
  if(missing(topology))
    stop("Missing parameter for confidence.fit.conditional function: confidence.fit.conditional(topology)", call. = FALSE)
  names <- rownames(topology@estimated.cond.probs)
  v <- as.vector(topology@estimated.cond.probs)
  names(v) <- names
  par("mar" = c((max(nchar(names))/1.6), 2, 2, 1), "mfrow" = c(1,1))
  barplot(v, las = 3, main = "Fit (conditional probability)")
  print(v)
}

#' @rdname confidence 
#' @usage confidence.single(topology)
#' @details \code{confidence.single} plot the difference between the observed and fitted  probability of each event
confidence.single <- function(topology){
  if(missing(topology))
    stop("Missing parameter for confidence.single function: confidence.single(topology)", call. = FALSE)
  data <- topology@marginal.probs
  fit <- topology@estimated.marginal.probs
  names <- rownames(data)
  v <- as.vector(data - fit)
  names(v) <- names
  par("mar" = c((max(nchar(names))/1.6), 2, 1, 1), "mfrow" = c(1,1))
  barplot(v, las = 3, main = "Confidence (events probability)")
  print(v)
}

#' @rdname confidence 
#' @usage confidence.joint(topology)
#' @details \code{confidence.joint} plot the pairwise difference between the observed and fitted joint probability of the events
confidence.joint <- function(topology){
  if(missing(topology))
    stop("Missing parameter for confidence.joint function: confidence.joint(topology)", call. = FALSE)
  data <- topology@joint.probs
  fit <- topology@estimated.joint.probs
  v <- data - fit
  print(v)
  levelplot(v, xlab = "", ylab = "", scales = list(x = list(alternating = 2, rot = 90), tck = 0), main = "Confidence (joint probability)")
}

#' @rdname confidence 
#' @usage confidence.conditional(topology)
#' @details \code{confidence.conditional} plot the pairwise difference between the observed and fitted conditional probability of the events
confidence.conditional <- function(topology){
  if(missing(topology))
    stop("Missing parameter for confidence.conditional function: confidence.conditional(topology)", call. = FALSE)
  data <- topology@cond.probs
  fit <- topology@estimated.cond.probs
  names <- rownames(data)
  v <- as.vector(data - fit)
  names(v) <- names
  par("mar" = c((max(nchar(names))/1.6), 2, 2, 1), "mfrow" = c(1,1))
  barplot(v, las = 3, main = "Confidence (conditional probability)")
  print(v)
}
