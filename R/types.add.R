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

#' @export types.add
#' @title add a new type of event (e.g., missense point mutation)
#'
#' @description
#' \code{types.add} sets a global data frame 'types' that contains all types defined.  Types can be added and refined incrementally, in any order.
#'
#' @details
#' \code{types.add} defines a type of event considered  at a time. If the type was previously defined, its definition is updated to keep track of its last definition. A consistency check is performed to ensure that the type is valid. Types must be defined before events are loaded.
#' 
#' @param type.name The type label. All type labels are strings.
#' @param color.name The type color. All R's color definitions are allowed.
#' @examples
#' types.add("gain", "red")
#' 
types.add <- function(type.name, color.name){
	
    if(missing(type.name) || missing(color.name)){
        stop("Missing parameter for types.add function: types.add(type.name, color.name)", call. = FALSE)
    }
    else{
      # If a global types variable is found, the new definition is queued to the definitions found.
      if(exists("types") && (length(types) > 0)){
      	types <- types
        types.local <- rbind(types, data.frame(type = type.name, color = color.name, stringsAsFactors = FALSE))
        
      }
      else{
      	types.local <- data.frame(type = type.name, color = color.name, stringsAsFactors = FALSE)
      }
     
      # The user is free to leave spaces between each element in the definition, definitions file is more clear this way.
      types.local <- trim.types(types.local)
    
      # The check function perform consistency and correctness checks.
      types.local <- check.types(types.local, FALSE)
      
      cat(paste("Set color value \"", color.name , "\" for events of type \"", type.name, "\"\n", sep =""))
      
      assign("types", types.local, envir = .GlobalEnv)
    }
    
 
}