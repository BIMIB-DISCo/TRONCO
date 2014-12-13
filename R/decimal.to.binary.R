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

#convert an integer decimal number to binary
#INPUT:
#num.decimal: decimal integer to be converted
#num.bits: number of bits to be used
#RETURN:
#num.binary: binary conversion of num.decimal
"decimal.to.binary" <-
function(num.decimal, num.bits) {
    #structure where to save the result
    num.binary = rep(0,num.bits);
    #convert the integer decimal number to binary
    pos = 0;
    while(num.decimal>0) {
        #compute the value of the current step
        num.binary[num.bits-pos] = num.decimal %% 2;
        #divide the number by 2 for the next iteration
        num.decimal = num.decimal %/% 2;
        pos = pos + 1;
    }
    return(num.binary);
}
