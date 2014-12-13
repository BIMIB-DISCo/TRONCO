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

#estimate the error rates by "L-BFGS-B" optimization in terms of L2-error
#INPUT:
#marginal.probs: marginal probabilities
#joint.probs: joint probabilities
#parents.pos: which event is the parent? 0 if none, a number otherwise
#RETURN:
#estimated.error.rates: estimated probabilities, false positive and false negative error rates
"estimate.tree.error.rates" <-
function(marginal.probs,joint.probs,parents.pos) {
    #function to be optimized by "L-BFGS-B" optimization in terms of L2-error
	f.estimation <- function(errors) {
        #set the current error rates with the starting point of the optimization being (e_pos,e_neg) = 0
        error.rates = list(error.fp=errors[1],error.fn=errors[2]);
        #estimate the observed probabilities given the error rates
        estimated.probs = estimate.tree.probs(marginal.probs,joint.probs,parents.pos,error.rates);
        #evaluate the goodness of the estimatione by L2-error on the estimated marginal and joint probabilities
        error.estimation = sum((marginal.probs-estimated.probs$marginal.probs)^2)+sum((joint.probs-estimated.probs$joint.probs)^2);
        return(error.estimation);
	}
    #the estimation is performed as in Byrd et al (1995)
    #this method allows for box constraints, i.e., each variable can be given a lower and/or upper bound
	estimated.error.rates = optim(c(0.00,0.00),f.estimation,method="L-BFGS-B",lower=c(0.00,0.00),upper=c(0.49,0.49))$par;
    #structure to save the results
    estimated.error.rates = list(error.fp=estimated.error.rates[1],error.fn=estimated.error.rates[2]);
    return(estimated.error.rates);
}
