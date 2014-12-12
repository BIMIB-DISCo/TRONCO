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

#estimate the probability of observing each sample in the dataset given the reconstructed topology
#INPUT:
#dataset: a valid dataset
#reconstructed.topology: the reconstructed topology
#estimated.marginal.probabilities: estimated marginal probabilities of the events
#estimated.conditional.probabilities: estimated conditional probabilities of the events
#error.rates: error rates for false positives and false negatives
#RETURN:
#probabilities: probability of each sample
"estimate.tree.samples" <- function(dataset, reconstructed.topology, estimated.marginal.probabilities, estimated.conditional.probabilities, error.rates) {
    #structure where to save the probabilities of the samples
    probabilities = array(-1,c(nrow(dataset),1));
    #topological properties:
    #1. tree number
    #2. parent
    #3. level in the tree
    topology.structure = array(0,c(nrow(reconstructed.topology),3));
    #go through the subtrees within the topology of four 
    tree.count = 0;
    for (i in 1:nrow(reconstructed.topology)) {
        #if node i has no parents, it is a root
        if(length(which(reconstructed.topology[,i]==1))==0) {
            tree.count = tree.count + 1;
            level = 1;
            #set the parameters for the root
            topology.structure[i,1] = tree.count;
            topology.structure[i,2] = -1;
            topology.structure[i,3] = level;
            curr.node = i;
            #go through this tree
            while (length(curr.node)>0) {
                #move to the next level
                level = level + 1;
                new.node = vector();
                for (j in 1:length(curr.node)) {
                    curr.new.node = which(reconstructed.topology[curr.node[j],]==1);
                    if(length(curr.new.node)>0) {
                        new.node = c(new.node,curr.new.node);
                        for (k in 1:length(curr.new.node)) {
                            #number of the current subtree
                            topology.structure[curr.new.node[k],1] = tree.count;
                            #parent of the current node
                            topology.structure[curr.new.node[k],2] = curr.node[j];
                            #level of this node
                            topology.structure[curr.new.node[k],3] = level;
                        }
                    }
                }
                curr.node = new.node;
            }
        }
    }
    #go through the dataset and evalutate the probability of each sample
    for (i in 1:nrow(dataset)) {
        sample.probability = 1;
        for (j in 1:tree.count) {
            #probability of this subtree (without any knowledge, I set it to 1)
            curr.sample.probability = 1;
            #entries referring to this subtree
            curr.entry = which(topology.structure[,1]==j);
            #samples of each element of this subtree
            curr.sample = dataset[i,curr.entry];
            #parents of each element of this subtree
            curr.parents = topology.structure[curr.entry,2];
            #level of each element of this subtree
            curr.levels = topology.structure[curr.entry,3];
            #set the probability as the one of the root of this tree
            curr.sample.probability = curr.sample.probability * estimated.marginal.probabilities[curr.entry[which(curr.levels==1,arr.ind=TRUE)],1];
            #set the maximum level of this subtree
            max.level = curr.levels[which.max(curr.levels)];
            #if I have at least one event in this sample
            if(length(curr.sample[curr.sample==1])>0) {
                #visit the nodes starting from the lower level
                is.valid = TRUE;
                for (k in max.level:1) {
                    curr.level.nodes = which(curr.levels==k,arr.ind=TRUE);
                    #if I'm not on a root
                    if(k>1) {
                        curr.level.samples = curr.sample[curr.level.nodes];
                        #if I have at least one event at this level
                        if(length(curr.level.samples[curr.level.samples==1])>0) {
                            #I can not have a child without its parent
                            curr.level.parent = curr.parents[curr.level.nodes];
                            for (p in 1:length(curr.level.parent)) {
                            	if(dataset[i,curr.level.parent[p]]==0 && dataset[i,curr.entry[curr.level.nodes[p]]]==1) {
                                	is.valid = FALSE;
                                	break;
                            	}
                            }
                        }
                        #if the sample is valid
                        if(is.valid==TRUE) {
                            #add the probability of each edge
                            curr.level.parent = curr.parents[curr.level.nodes];
                            for (p in 1:length(curr.level.parent)) {
                            	if(dataset[i,curr.level.parent[p]]==1 && dataset[i,curr.entry[curr.level.nodes[p]]]==0) {
                                    curr.sample.probability = curr.sample.probability * (1 - estimated.conditional.probabilities[curr.entry[curr.level.nodes[p]],1]);
                            	}
                                else if(dataset[i,curr.level.parent[p]]==1 && dataset[i,curr.entry[curr.level.nodes[p]]]==1) {
                                    curr.sample.probability = curr.sample.probability * estimated.conditional.probabilities[curr.entry[curr.level.nodes[p]],1];
                                }
                            }
                        }
                    }
                    if(is.valid==FALSE) {
                        curr.sample.probability = 0;
                        break;
                    }
                }
                if(is.valid==FALSE) {
                    sample.probability = 0;
                    break;
                }
            }
            #if this sample has no events for this tree
            else {
                curr.sample.probability = 1 - curr.sample.probability;
            }
            #update the probability of the topology with the one of this sample
            sample.probability = sample.probability * curr.sample.probability;
            if(sample.probability==0) {
                break;
            }
        }
        probabilities[i,1] = sample.probability;
    }
    #correct the estimation by the error rates
    errors.matrix <- array(0,c(nrow(probabilities),nrow(dataset)));
    for (i in 1:nrow(probabilities)) {
        for (j in 1:nrow(dataset)) {
            curr.sample.x = as.numeric(dataset[i,]);
            curr.sample.y = as.numeric(dataset[j,]);
            errors.matrix[i,j] = (1-error.rates$error.fp)^((1-curr.sample.x)%*%(1-curr.sample.y))*error.rates$error.fp^((1-curr.sample.x)%*%curr.sample.y)*(1-error.rates$error.fn)^(curr.sample.x%*%curr.sample.y)*error.rates$error.fn^(curr.sample.x%*%(1-curr.sample.y));
        }
    }
    probabilities[,1] = as.numeric(as.vector(probabilities)%*%errors.matrix);
    return(probabilities);
}
