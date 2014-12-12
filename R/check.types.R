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

# Checks if types variable is well declared 
check.types <- function(types, file.or.function = FALSE){
  
  if(file.or.function)
    msg <- "file"
  else
    msg <- ""
    
  # Looks for declaration duplicates
  dup <- c()
  # K flag is set if duplicate by type name are found 
  k <- FALSE
  # P plag is set if duplicates color are found
  p <- FALSE
  
  
  i <- 1
  while(i <= nrow(types)){
    j <- 1
    while(j <= (i-1)){
      # Checks if a type is duplicated by type name
      if((types[j, "type"] == types[i, "type"])){
        dup <- union(dup, c(j))
        k <- TRUE
      }
      j <- j + 1
    }
    i <- i + 1
  }
 
  # If any duplicate is found the lastone is kept and a warning message is shown
  if(k){
  	types.dup <- types
  	types <- types[-1*dup,]
    shown <- c("")
    for(i in 1:length(dup)){
      if(! any(shown == types.dup[dup[i], "type"])){
      	r <- match(types.dup[dup[i], "type"], types[,"type"])
        warning(paste("Event type ",
        			  toString(types.dup[dup[i], "type"]),
        			  " redefined, now has color: ", 
        			  toString(types[r,"color"]), 
        			  sep = ""), call. = FALSE)
        			  
        shown <- union(shown, types.dup[dup[i], "type"])
      }
    }
  }
  
  i <- 1
  while(i <= nrow(types)){
    j <- 1
    while(j <= (i-1)){
      # Checks if any set of type have the same color
      if((types[j, "color"] == types[i, "color"]) && (types[j, "type"] != types[i, "type"])){
        p <- TRUE
      }
      j <- j + 1
    }
    i <- i + 1
  }
  # If there is more than one type with the same color declared a warning message is displayed 
  if(p)
    warning("There are multiple events with the same color defined.", call. = FALSE)
  
  
  # reset
  t <- types
  row.names(t) <- c(1:nrow(t))
  return(t)
   
}