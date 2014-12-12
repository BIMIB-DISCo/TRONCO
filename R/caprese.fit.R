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

#reconstruct the best tree-like topology
#INPUT:
#dataset: a dataset describing a progressive phenomenon
#lambda: shrinkage parameter (value in [0,1])
#verbose: should I print the warnings? Yes if TRUE, no otherwise
#RETURN:
#topology: the reconstructed tree-like topology
"caprese.fit" <-
function(dataset,lambda,verbose) {
	#structure to compute the observed marginal and joint probabilities
	pair.count <- array(0, dim=c(ncol(dataset), ncol(dataset)));
	#compute the probabilities on the dataset
 	for(i in 1:ncol(dataset)) {
		for(j in 1:ncol(dataset)) {
        	val1 = dataset[ ,i];
        	val2 = dataset[ ,j];
            pair.count[i,j] = (t(val1) %*% val2);
        }
	}
    #marginal.probs is an array of the observed marginal probabilities
    marginal.probs <- array(as.matrix(diag(pair.count)/nrow(dataset)),dim=c(ncol(dataset),1));
    #joint.probs is an array of the observed joint probabilities
    joint.probs <- as.matrix(pair.count/nrow(dataset));
    #reconstruct the causal topology
    best.parents = get.tree.parents(marginal.probs,joint.probs,lambda);
    #create the structures where to save the results
    parents.pos <- best.parents$parents;
    conditional.probs <- array(-1, dim=c(length(parents.pos),1));
    adj.matrix <- array(0, dim=c(length(parents.pos),length(parents.pos)));
    #set the parents names and the structures
    for(i in 1:ncol(dataset)) {
        #if the node has a parent
        if(parents.pos[i,1]!=-1) {
            #Note: [i,j] = 1 means that i is causing j
            adj.matrix[parents.pos[i,1],i] = 1;
            #compute the conditional probability of P(CHILD=1|PARENT=1)
        	conditional.probs[i,1] = best.parents$joint.probs[parents.pos[i,1],i]/best.parents$marginal.probs[parents.pos[i]];
        }
        #if the node has no parent, its conditional probability is set to 
        else {
        	conditional.probs[i,1] = 1;
    	}
	}

	#estimate the error rates and, given them, the probabilities
    estimated.error.rates = estimate.tree.error.rates(best.parents$marginal.probs,best.parents$joint.probs,parents.pos);
    estimated.probabilities = estimate.tree.probs(best.parents$marginal.probs,best.parents$joint.probs,parents.pos,estimated.error.rates);
    #structures where to save the probabilities
    probabilities = list(marginal.probs=best.parents$marginal.probs,joint.probs=best.parents$joint.probs,conditional.probs=conditional.probs,estimated.marginal.probs=estimated.probabilities$marginal.probs,estimated.joint.probs=estimated.probabilities$joint.probs,estimated.conditional.probs=estimated.probabilities$conditional.probs);
    parameters = list(lambda=lambda);
    #return the results
    topology = list(dataset=dataset, probabilities=probabilities,error.rates=estimated.error.rates,pr.score=best.parents$pr.score,adj.matrix=adj.matrix,parameters=parameters);
    return(topology);
}
