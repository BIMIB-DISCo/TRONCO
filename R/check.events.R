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

# Checks if events variable is well declared 
check.events <- function(events, types, file.or.function = FALSE){
 
  # if file.or.function is set TRUE the error messages are about
  # the function that loads events from file 
  if(file.or.function)
    msg <- "file"
  else
    msg <- ""

  # Looks for declaration duplicates
  dup <- c()
  # p flag is set if duplicate by name or type are found
  p <- FALSE
  i <- 1
  while(i <= nrow(events)){
    j <- 1
    while(j <= (i-1)){
      # Checks if an event is duplicated by name and type
      if((events[j, "event"] == events[i, "event"]) && (events[j, "type"] == events[i, "type"])){          
         dup <- union(dup, c(j))
         p <- TRUE
      }
      j <- j + 1
    }
    i <- i + 1
  }
  if(p){
  	# Displays once duplicate messages
    shown <- c("")
    for(i in 1:length(dup)){
      if(! any(shown == events[dup[i], "event"])){
        warning(paste(toString(events[dup[i], "event"])," event redefined!", sep = ""), call. = FALSE)
        shown <- union(shown, events[dup[i], "event"])
      }
    }
  }
  
  # k flag is set if more than one event is assigned to one column
  dup.col <- c()
  k <- FALSE
  i <- 1
  while(i <= nrow(events)){
    j <- 1
    while(j <= (i-1)){
      if(events[j, "column"] == events[i, "column"]){
        dup.col <- union(dup.col, c(j))
        # If k flag is set warning messages is shown and the event will be discarded
        k <- TRUE
      }
      j <- j + 1
    }
    i <- i + 1
  }
  if(k){
    for(i in 1:length(dup.col))
      #warning(paste("Event ", toString(events[dup[i], "event"]), " of type ", toString(events[dup[i], "type"]), " is redefined", sep = ""), call. = FALSE)
      warning(paste("Event ", toString(events[dup.col[i], "event"]), " of type ", toString(events[dup.col[i], "type"]), " is in a column redefined below, row discarded", sep = ""), call. = FALSE)
    
 }
 if(p || k){
 	dup <- union(dup, dup.col)
 	events <- events[-1*dup,]
 }
 
 
 # If the types variable is found as it must, for each row this section checks the integrity of each column
 # in the events variable
 if(exists("types") && (length(types) > 0)){
    for(i in 1:nrow(events)){
      if(!search.type(toString(events[i, "type"])))
        warning(paste("Type ", toString(events[i,"type"]), " is not found!", sep = ""), call. = FALSE)
      if(!is.wholenumber(events[i, "column"]))
        stop(paste(toString(events[i,"event"])," column must be an integer value!", sep = ""), call. = FALSE)
      if(events[i, "column"] < 0){
        e <- toString(events[i,"event"])
        events <- events[-1*i,]
        stop(paste(e," column must be an integer positive value row discarded!", sep = ""), call. = FALSE)	
      }
    }
  }else
    stop("types variable not found!", call. = FALSE)
  
  e <- events
  row.names(e) <- c(1:nrow(e))
  return(e)
  
}