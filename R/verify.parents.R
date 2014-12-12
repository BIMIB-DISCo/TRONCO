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

#verify the independent progression filter
#INPUT:
#best.parents: best edges to be verified
#marginal.probs: observed marginal probabilities
#joint.probs: observed joint probabilities
#RETURN:
#best.parents: list of the best valid parents
"verify.parents" <-
function(best.parents,marginal.probs,joint.probs) {
	#verify the condition for the best parent of each node
	for (i in 1:length(best.parents)) {
        #if there is a connection, i.e. the node is not already attached to the Root
		if(best.parents[i]!=-1) {
            #score for the root as the parent of this node
			w.root.node = 1/(1+marginal.probs[i]);
			#compute the scores for the edges to all the other upstream nodes
			attach.to.root = 1;
			for (j in 1:length(marginal.probs)) {
                #if the connection is valid and the parent node has greater probability
                #i.e. it is before the child in temporal order
				if(i!=j && marginal.probs[j]>marginal.probs[i]) {
					w.parent.node = (marginal.probs[j]/(marginal.probs[i]+marginal.probs[j]))*(joint.probs[i,j]/(marginal.probs[i]*marginal.probs[j]));
                    #the parent is valid if this condition is valid at least one time (i.e. for at least one of the upstream nodes)
                    #meaning that if we find out that a connection is not spurious for any node, the best parent is not spurious as well
					if(w.root.node<=w.parent.node) {
						attach.to.root = 0;
						break;
					}
				}
			}
			#connect the node to the Root if the flag is true
			if(attach.to.root==1) {
				best.parents[i] = -1;
			}
		}
	}
	return(best.parents);
}
