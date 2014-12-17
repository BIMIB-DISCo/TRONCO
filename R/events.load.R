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

#' @export events.load
#' @title load a set of events from file
#'
#' @description
#' \code{events.load} sets a global data frame 'events' that contains all event definitions found in a specified file or dataset to be validated. This is a way to automatise calls to function \code{events.add} for a bunch of events.
#'
#' @details
#' \code{events.load} load a set of events from a given file. The input file must be structured as a CSV file, where each event is defined on a separate line in the format: eventName, typeName, columnNumber.
#' 
#' @seealso \code{\link{events.add}}
#' @param data.input The input file path or a dataset to be validated.
#' 
events.load <- function(data.input){
	
  err <- ""
  message <- "The definition file contains errors!"
  
  if(missing(data.input))
    stop("Missing parameter for events.load function: events.load(data.input)", call. = FALSE)
  
  	 # Types must be defined before the event definition
    if(exists("types") && (length(types) > 0)){
    	
      types <- types
      # If a global events variable is found, the new definition is queued to the definitions found.
      if(exists("events") && (length(events) > 0)){
      	events <- events
        
      	if(is.data.frame(data.input))
      	  events.file <- data.input
        else{
          # If the pathname is correct
          if(file.exists(data.input)){
          	# Definition file may contain error such as the lack of a comma or columns, a try-catch manages this.
            err <- tryCatch( events.file <- suppressWarnings(read.table(data.input, sep = ",", col.names = c("event", "type", "column"), stringsAsFactors = FALSE)),
                      error = function(e) err <- message)
            if(toString(err) == message)
              stop(err, call. = FALSE)
          }else
            stop("File not found!", call. = FALSE)
        }
      	if(nrow(events.file) > 0)
      	  # The new definitions are queued.   
      	  events.local <- rbind(events, events.file)
        else
      	  stop("Empty set of events at input file path or dataset!", call. = FALSE)
      }
      else{
        if(is.data.frame(data.input))
          events.local <- data.input
        else{
          # If the pathname is correct
          if(file.exists(data.input)){
            err <- tryCatch( events.local <- suppressWarnings(read.table(data.input, sep = ",", col.names = c("event", "type", "column"), stringsAsFactors = FALSE)),
                      error = function(e) err <- message)
            if(toString(err) == message)
              stop(err, call. = FALSE)
          }else
            stop("File not found!", call. = FALSE)
        }
        if(nrow(events.local) == 0)
          stop("Empty set of events at input file pathor dataset!", call. = FALSE)
      }
      
      # The user is free to leave spaces between each element in the definition, definitions file is more clear this way.
      events.local <- trim.events(events.local)
      
      
      # The check function perform consistency and correctness checks.
      events.local <- check.events(events.local, types, TRUE)
      
      
      for(i in 1:nrow(events.local))
        cat(paste("Added event \"", events.local[i, "event"] , "\" of type \"", events.local[i, "type"], "\" (color: \"", 
        	toString(search.type.info(events.local[i, "type"])[,"color"]),"\"), dataset column \"", events.local[i, "column"],"\"\n", sep =""))
      
      assign("events", events.local, envir = .GlobalEnv)

    }else
      stop("Types not defined!", call. = FALSE)
  
}