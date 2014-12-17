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

#' @export events.add
#' @title add a new event (e.g., a missense point mutation for EGFR)
#'
#' @description
#' \code{events.add} sets a global data frame 'events' that contains all the events defined. Events can be added and refined incrementally, in any order.
#'
#' @details
#' \code{events.add} allows to define one event at a time. If the event was previously defined, its definition is updated to keep track of its last definition. A consistency check is performed to ensure that the type of defined event is valid. Thus, types must be defined before events are loaded (see \code{types.add}, \code{types.load}).
#' 
#' @param event.name The event label(e.g., 'EGFR') . All event labels are strings.
#' @param type.name The type name of this event (e.g., 'missense point'). Type names must refer to types loaded before adding an event, a consistency check raises an error if the type name is unknown.
#' @param column.number The dataset column to which this event is associated. Column number must be an
#' integer positive value.
#' @examples
#' types.add("gain", "red")
#' events.add("8q+", "gain", 1)
#' 
events.add <- function(event.name, type.name, column.number = NA){
  
  if(missing(event.name) || missing(type.name) || missing(column.number))
    stop("Missing parameter for events.add function: events.add(event.name, type.name, column.number)")
  
  # Types must be defined before the event definition
  if(exists("types") && (length(types) > 0)){
  	types <- types
  	
  	# Column number has to be a valid number
    if(!is.na(column.number)){
    	
    	# Column number has to be an integer number
    	if(is.wholenumber(column.number)){
      		
      		# Each number given in input from the console is considered to be a "numeric" value.
      		c <- as.integer(column.number)
      		
      		# If a global events variable is found, the new definition is queued to the definitions found.
      		if(exists("events") && (length(events) > 0)){
      			events <- events
        		events.local <- rbind(events, data.frame(event = event.name, type = type.name, column = c, stringsAsFactors = FALSE))
        	}
      		else
        		events.local <- data.frame(event = event.name, type = type.name, column = c, stringsAsFactors = FALSE)
          
    
          	# The user is free to leave spaces between each element in the definition, definitions file is more clear this way.
      		events.local <- trim.events(events.local)
      		
      		# The check function perform consistency and correctness checks.
      		events.local <- check.events(events.local, types, FALSE)
          
          	cat(paste("Added event \"", event.name , "\" of type \"", type.name, "\" (color: \"", 
          	  	toString(search.type.info(type.name)[,"color"]),"\"), dataset column \"", column.number,"\"\n", sep =""))
          	  	
          	assign("events", events.local, envir = .GlobalEnv)
          	  
    	}
      else
    	  stop("Column must be an integer value!", call. = FALSE)
    }
    else
      stop("Cannot create an event without assigning it a column number.", call. = FALSE)  
  }
  else
    stop("types  variable not defined!", call. = FALSE)
}