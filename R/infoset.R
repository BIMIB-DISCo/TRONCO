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

infoset <- function(merged.columns, deleted.column, verbose = FALSE){
  
  events <- events
  # Gets all useful info for each node and returns a matrix that 
  # contains color and label to be assigned to each event.
  get.all.info  <- function(events, types){
    all.vis.labels <- c()
    all.labels  <- c()
    colors <- c()
    for(i in 1:nrow(events)){
      
      # For the event in column i look for color, label name for computation
      # and label good to graph drawing 
      color  <- search.type.info(search.event(i)$type)$color
      label <- paste(toString(search.event(i)$event), toString(search.event(i)$type), sep=":")
      label.vis <- toString(search.event(i)$event)
		
	  if(is.null(color))
	  	stop(paste("For the event \"", search.event(i)$event, "\" no color can be associated to the given type, check the types definition!", sep = ""), call. = FALSE)
      # Load
      colors <- c(colors, color)
      all.vis.labels  <- c(all.vis.labels, label.vis)
      all.labels <- c(all.labels, label)
    }
    
    all.info <- list()
    all.info$all.labels  <- all.labels
    all.info$all.vis.labels <- all.vis.labels
    all.info$colors <- colors
    return(all.info)
    
  }
  
  all.info <- get.all.info(events, types)
  
  # Events are merged as merged.columns indicates.
  if(length(merged.columns) > 0 && !is.na(merged.columns)){
  	
    for(i in 1:nrow(merged.columns)){
    
      merged <- merged.columns[i,1]
      with <- merged.columns[i,2]
      
      # all.labels contains all the "computations label", which are composed by event_name:type, 
      # this turns the label into a key for the computation process
      all.info$all.labels[merged] <- paste(all.info$all.labels[merged], all.info$all.labels[with], sep = " - ")
      
      # all.vis.labels cotains all labels to be assigned to each node in the final plot
      all.info$all.vis.labels[merged] <- paste(all.info$all.vis.labels[merged], all.info$all.vis.labels[with], sep = " - ")

      if(verbose)
        cat(paste("Column ", toString(merged), " event: ", search.event(merged)$event, " merged with column ", 
        	toString(with), " event: ", search.event(with)$event, " ", search.event(deleted.column[i])$eventcolor, 
        	"lightyellow color", " is assigned\n", sep = ""))
        	
      all.info$colors[merged] <- "lightyellow"
      
      
    }
    if(! search.type("merged"))
      types <- rbind(types, data.frame(type = "merged", color = "lightyellow", stringsAsFactors = FALSE))
      assign("types", types, envir = .GlobalEnv)
  }
  
  # Deletes columns if some column number is in deleted.column vector.
  # WARINIG: All columns merged and useless are already put in the deleted.column vetor by CAPRESE algorithm
  if(length(deleted.column) > 0 && !is.na(deleted.column[1])){
    all.info$all.labels <- all.info$all.labels[-1*deleted.column]
    all.info$all.vis.labels <- all.info$all.vis.labels[-1*deleted.column]
    all.info$colors <- all.info$colors[-1*deleted.column]
    if(verbose)
      for(i in 1:length(deleted.column))
        cat(paste("Column ", toString(deleted.column[i]), " type: ", 
                  search.event(deleted.column[i])$type, " event: ",
                  search.event(deleted.column[i])$event,
                  " deleted\n", sep = ""))
  }
  
  return(all.info)
  
}