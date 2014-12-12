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

#compute the probability raising based scores
#INPUT:
#marginal.probs: observed marginal probabilities
#joint.probs: observed joint probabilities
#lambda: shrinkage parameter (value between 0 and 1)
#RETURN:
#scores: probability raising based scores
"get.tree.scores" <-
function(marginal.probs,joint.probs,lambda) {
    #structure where to save the probability raising scores
	pr.score = array(-1, dim=c(nrow(marginal.probs),nrow(marginal.probs)));
	#compute the probability raising based scores
	for (i in 1:ncol(pr.score)) {
        for (j in 1:ncol(pr.score)) {
            #alpha is the probability raising model of causation (raw model estimate)
            alpha = ((joint.probs[i,j]/marginal.probs[i])-((marginal.probs[j]-joint.probs[i,j])/(1-marginal.probs[i])))/((joint.probs[i,j]/marginal.probs[i])+((marginal.probs[j]-joint.probs[i,j])/(1-marginal.probs[i])));
            #beta is the correction factor (based on time distance in terms of statistical dependence)
            beta = (joint.probs[i,j]-marginal.probs[i]*marginal.probs[j])/(joint.probs[i,j]+marginal.probs[i]*marginal.probs[j]);
            #the overall estimator is a shrinkage-like combination of alpha and beta
            #the scores are saved in the convention used for an ajacency matrix, i.e. [i,j] means causal edge i-->j
            pr.score[i,j] = (1-lambda)*alpha + lambda*beta;
        }
	}
	scores = list(marginal.probs=marginal.probs,joint.probs=joint.probs,pr.score=pr.score);
	return(scores);
}
