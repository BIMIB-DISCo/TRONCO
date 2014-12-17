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

#' @export reset.events
#' @export reset.types
#' @export reset
#' @name reset
#' @title reset
#' @description
#' A set of functions to reset events, types and data.values variables
#' 
#' @usage reset.events() 
#' @details \code{reset.events} Resets the events variable
#' @examples
#' reset.events()
reset.events <- function(){
   assign("events", NULL, envir = .GlobalEnv)
}
#' @rdname reset
#' @usage reset.types()
#' @details \code{reset.types()} Resets the types variable
#' @examples
#' reset.types()
reset.types <- function(){
   assign("types", NULL, envir = .GlobalEnv)
}
reset.data.values <- function(){
	assign("data.values", NULL, envir = .GlobalEnv)
}
#' @rdname reset
#' @usage reset()
#' @details \code{reset()} Resets types, events and data.values variables
#' @examples
#' reset()
reset <- function(){
  reset.events()
  reset.types()
  reset.data.values()
  if(exists("num.hypotheses")) {
  	assign("num.hypotheses", 0, envir = .GlobalEnv)
  }
  if(exists("llist")) {
  	assign("llist", vector(), envir = .GlobalEnv)
  }
  if(exists("hlist")) {
  	assign("hlist", vector(), envir = .GlobalEnv)
  }
}
