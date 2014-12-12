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

#' @export data.load
#' @title load a dataset (binary matrix) from a file or a preloaded dataset.
#'
#' @description
#' \code{data.load} sets a global data frame 'data.values' that contains the dataset loaded from an input file.
#'
#' @details
#' \code{data.load} loads a dataset from disk and associates all columns in the dataset to a specified event. Thus, types and events must be specified before calling this function to ensure a consistency check is performed on the input dataset (see \code{types.load}, \code{types.add},  \code{events.load}, \code{events.add} to load/add types/events).
#' 
#' @param data.input The input file path. or a dataset loaded by \code{data} function
#' 
 
data.load <- function(data.input){
  
  if(missing(data.input))
    stop("Missing parameter for data.load function: data.load(data.input)", call. = FALSE)
  # To load data from the file linked by the given path the variables events and types must be already defined  
  if(exists("types") && exists("events") && (length(types) > 0) && (length(events) > 0)){
    events <- events
  	types <- types
  	
  	# The read.table function loads data in the data.values variable
  	# WARNING: warnings for this instruction are suppressed because they occurs even when 
  	# the user do not leave a final return in the definitions file, we doesn't matter this.
    if(!is.data.frame(data.input))
    	data.values <- suppressWarnings(read.table(data.input, header = FALSE))
    else
    	data.values <- data.input

    # The number of elements found in events variable must be enough to make an association with the dataset columns.
    if(ncol(data.values) > nrow(events)) stop("Events definition is not complete, not enough events defined!")
    
    # If an event declaration contains a number of column out of bound, an error message is displayed.
    for(i in 1:nrow(events))
      if(events[i, "column"] > ncol(data.values))
        stop(paste("Found an event with column number greater than dataset colum number, event:", events[i, "event"]))
    
    column.names <- c()
    
    # Each event is associated to his column.
    for(i in 1:ncol(data.values)){
      event <- search.event(i)
      column.name <- c(paste(toString(event$event), " (", toString(event$type), ", column ", toString(event$column), ")", sep = ""))
      column.names <- c(column.names, column.name)
    }
   
    # check datasets
	  # colnames(data.values) <- column.names
	  tmpdata <- check.dataset(data.values, FALSE)
    # Collects info for the object of class topology.
    info  <- infoset(tmpdata$invalid.events$merged.events, tmpdata$invalid.events$removed.events, verbose = TRUE)
    assign("invalid.events", tmpdata$invalid.events , envir = .GlobalEnv)
    data.values <- tmpdata$dataset
    
    num.events <- ncol(data.values)
    num.samples <- nrow(data.values)
    num.hypotheses <- 0
    
    assign("num.events", num.events, envir = .GlobalEnv)
    assign("num.samples", num.samples, envir = .GlobalEnv)
    assign("num.hypotheses", num.hypotheses, envir = .GlobalEnv)
    
    if(!is.null(nrow(tmpdata$invalid.events$merged.events))){
    	cat("\nThe available events are:\n")
    	events.label <- info$all.vis.labels
    	out.matrix <- events
      
    	merged.type <- tmpdata$invalid.events$merged.events[,1]
    	out.matrix[merged.type,"type"] <- "merged"
      
      out.matrix <- out.matrix[-1*tmpdata$invalid.events$removed.events,]
      
      names <- cbind(info$all.vis.labels)
      out.matrix[,"event"] <- names
      
      #columns numbers
      out.matrix[, "column"] <- cbind(1:nrow(out.matrix))
    	row.names(out.matrix) <- c(1:nrow(out.matrix))
      
      print(out.matrix)
    		
    }
    
	colnames(data.values) <- info$all.labels
	assign("data.values", data.values, envir = .GlobalEnv)
	if(!is.data.frame(data.input))
    cat(paste("Data successfully loaded and validated, dataset: ", toString(data.input),"\n", sep =""))
  else
    cat("Data frame validated and data.values variable is loaded in global environment")
    
  }
  else
    stop("Events or types variable not found: Complete the definitions!")
  
}
