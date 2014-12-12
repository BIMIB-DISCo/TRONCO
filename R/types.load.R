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

#' @export types.load
#' @title load a set of types from file
#'
#' @description
#' \code{types.load} sets a global data frame 'types' that contains all type definitions found in a specified file or dataset
#' to be validated.
#'
#' @details
#' \code{types.load} allows to load type definitions from a given file path. The file which contains
#' all the definitions must be structured as a csv file. All definitions are couple of values 
#' type name and color name as shown below:
#' 
#' typeName, colorName
#' ...     , ...
#' 
#' @seealso \code{\link{types.add}}
#' @param data.input The input file path or a dataset to be validated.

types.load <- function(data.input){
	
  err <- ""
  message <- "The definition file contains errors!"
  
  if(missing(data.input)){
    stop("Missing parameter for types.load function: types.load(data.input)", call. = FALSE)
  }
  
  
  	
  	# If a global events variable is found, the new definition is queued to the definitions found.
    if(exists("types") && (length(types) > 0)){
    	types <- types
      
    	if(is.data.frame(data.input))
    	  types.file <- data.input
    	else{
    	  # If the pathname is correct
    	  if(file.exists(data.input)){
          # Definition file may contain error such as the lack of a comma or columns, a try-catch manages this.
          err <- tryCatch( types.file <- suppressWarnings(read.table(data.input, sep = ",", col.names = c("type", "color"), stringsAsFactors = FALSE)),
                    error = function(e) err <- message)
          if(toString(err) == message)
            stop(err, call. = FALSE)
    	  }else
    	    stop("File not found!", call. = FALSE) 
    	}
      if(nrow(types.file) > 0)
         # The new definitions are queued.
        types.local <- rbind(types, types.file)
      else
    	  stop("Empty set of types at input file path or dataset!", call. = FALSE)
    }
    else{
      if(is.data.frame(data.input))
        types.local <- data.input
      else{
        # If the pathname is correct
        if(file.exists(data.input)){
          err  <- tryCatch( types.local <- suppressWarnings(read.table(data.input, sep = ",", col.names = c("type", "color"), stringsAsFactors = FALSE)),
                    error = function(e) err <- message)
          if(toString(err) == message)
            stop(err, call. = FALSE)
        }else
          stop("File not found!", call. = FALSE) 
      }
      if(nrow(types.local) == 0)
        stop("Empty set of types at input file path or dataset!", call. = FALSE)
    }

    # The user is free to leave spaces between each element in the definition, definitions file is more clear this way.
    types.local <- trim.types(types.local)
  
    # The check function perform consistency and correctness checks.
    types.local <- check.types(types.local, TRUE)
    
    for(i in 1:nrow(types.local))
      cat(paste("Set color value \"", types.local[i,"color"] , "\" for events of type \"", types.local[i,"type"], "\"\n", sep =""))
    
    assign("types", types.local, envir = .GlobalEnv)
   
}