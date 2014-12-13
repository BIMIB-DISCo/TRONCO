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

# Check if there is a type with the given tipe.name in type variable.
search.type <- function(type.name){
	types <- types
  	for(i in 1:nrow(types)){
    	if(toString(types[i,"type"]) == type.name)
    	  return(TRUE)
 	 }
  	return(FALSE)
}

# Search type information for the given type.name.
search.type.info <- function(name.type){
	types <- types
  	for(j in 1:nrow(types))
    	if(types[j,]$type == name.type)
    		return(types[j,])
}

# Search event informations for the given column.
search.event <- function(column.index){
	events <- events
  	column.index <- as.integer(column.index)
  	for(j in 1:nrow(events))
    	if(events[j,]$column == column.index)
      		return(events[j,])
  	stop(paste("Events definition is not complete!: events for column ", toString(column.index), " not found!", sep = ""), call.= FALSE)
}

exists.event <- function(event.name){
	events <- events
  	for(j in 1:nrow(events))
    	if(events[j,]$event == event.name)
      		return(events[j,])
  	return(NULL)
}

