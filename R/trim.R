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

# Resizes a string if spaces in the beginning and in the end are found
trim <- function(string){
  
  dim <- as.integer(nchar(string))
  str <- substr(string,1,1)
  p <- " "
  
  # Cuts the left blank space if at least one is found.
  while((str == p)){
    string <- substr(string,2,dim)
    dim <- as.integer(nchar(string))
    str <- substr(string,1,1)
  }
  
  dim <- as.integer(nchar(string))
  str <- substr(string,dim,dim)
  p <- " "
  
  # Cuts the right blank space if at least one is found.
  while((str == p)){
    string <- substr(string,1,(dim - 1))
    dim <- as.integer(nchar(string))
    str <- substr(string,dim,dim)
  }
  
  return(string) 
}