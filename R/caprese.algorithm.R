##################################################################################
#                                                                                #
# TRONCO: a tool for TRanslational ONCOlogy                                      #
#                                                                                #
##################################################################################
# Copyright (c) 2015, Marco Antoniotti, Giulio Caravagna, Luca De Sano,          #
# Alex Graudenzi, Ilya Korsunsky, Mattia Longoni, Loes Olde Loohuis,             #
# Giancarlo Mauri, Bud Mishra and Daniele Ramazzotti.                            #
#                                                                                #
# All rights reserved. This program and the accompanying materials               #
# are made available under the terms of the GNU GPL v3.0                         #
# which accompanies this distribution                                            #
#                                                                                #
##################################################################################


#' reconstruct the best tree-like topology
#' @title caprese.fit
#' @param dataset a dataset describing a progressive phenomenon
#' @param lambda shrinkage parameter (value in [0,1])
#' @param do.estimation should I perform the estimation of the error rates and probabilities?
#' @param silent execute the algorithm in silent mode
#' @return topology: the reconstructed tree-like topology
caprese.fit = function( dataset, lambda = 0.5 , do.estimation = FALSE, silent = FALSE ) {
	
	# start the clock to measure the execution time
	ptm <- proc.time();
	
	# structure with the set of valid edges
    # I start from the complete graph, i.e., I have no prior and all the connections are possibly causal
    adj.matrix = array(1,c(ncol(dataset),ncol(dataset)));
    colnames(adj.matrix) = colnames(dataset);
    rownames(adj.matrix) = colnames(dataset);
    
    # the diagonal of the adjacency matrix should not be considered, i.e., no self cause is allowed
    diag(adj.matrix) = 0;
    
    # check if the dataset is valid
    valid.dataset = check.dataset(dataset,adj.matrix,FALSE);
    adj.matrix = valid.dataset$adj.matrix;
	
	# marginal.probs is an array of the observed marginal probabilities
    marginal.probs <- valid.dataset$marginal.probs;
    # joint.probs is an array of the observed joint probabilities
    joint.probs <- valid.dataset$joint.probs;
    
    # reconstruct the causal topology
    best.parents = get.tree.parents(adj.matrix,marginal.probs,joint.probs,lambda);
    
    # create the structures where to save the results
    parents.pos <- best.parents$parents;
    conditional.probs <- array(-1, dim=c(length(parents.pos),1));
    adj.matrix <- array(0, dim=c(length(parents.pos),length(parents.pos)));
    colnames(adj.matrix) = colnames(dataset);
    rownames(adj.matrix) = colnames(dataset);
    confidence <- array(list(), c(3,1));
    confidence[[1,1]] = array(0, dim=c(length(parents.pos),length(parents.pos)));
    confidence[[2,1]] = best.parents$pr.score;
    confidence[[3,1]] = array(0, dim=c(length(parents.pos),length(parents.pos)));
    
    # set the parents names and the structures
    hypergeometric.pvalues = vector();
    for(i in 1:ncol(dataset)) {
        # if the node has a parent
        if(parents.pos[i,1]!=-1) {
            # Note: [i,j] = 1 means that i is causing j
            adj.matrix[parents.pos[i,1],i] = 1;
            # compute the conditional probability of P(CHILD=1|PARENT=1)
            conditional.probs[i,1] = best.parents$joint.probs[parents.pos[i,1],i]/best.parents$marginal.probs[parents.pos[i]];
        }
        # if the node has no parent, its conditional probability is set to 
        else {
			conditional.probs[i,1] = 1;
        }
        # compute the hypergeometric test
        for (j in i:ncol(dataset)) {
        		if(i!=j) {
        			
        			# confidence in temporal priority
        			confidence[[1,1]][i,j] = min(best.parents$marginal.probs[i]/best.parents$marginal.probs[j],1);
        		confidence[[1,1]][j,i] = min(best.parents$marginal.probs[j]/best.parents$marginal.probs[i],1);
        			
        			# compute the confidence by hypergeometric test
        			confidence[[3,1]][i,j] = phyper(best.parents$joint.probs[i,j]*nrow(dataset),best.parents$marginal.probs[i]*nrow(dataset),nrow(dataset)-best.parents$marginal.probs[i]*nrow(dataset),best.parents$marginal.probs[j]*nrow(dataset),lower.tail=FALSE);
        			confidence[[3,1]][j,i] = confidence[[3,1]][i,j];
        			# save all the valid pvalues
        			hypergeometric.pvalues = append(hypergeometric.pvalues,confidence[[2,1]][i,j]);
        			
        		}
        }
    }
	
    if(do.estimation) {
		#estimate the error rates and, given them, the probabilities
		estimated.error.rates = estimate.tree.error.rates(best.parents$marginal.probs,best.parents$joint.probs,parents.pos);
		estimated.probabilities = estimate.tree.probs(best.parents$marginal.probs,best.parents$joint.probs,parents.pos,estimated.error.rates);
	}
	else {
		estimated.error.rates = list(error.fp=NA,error.fn=NA);
		estimated.probabilities = list(marginal.probs=NA,joint.probs=NA,conditional.probs=NA);
	}
    error.rates = estimated.error.rates;
    estimated.probabilities.fit = estimated.probabilities;
	
    #structures where to save the results
    model = list();
    adj.matrix.fit = list();
    adj.matrix.fit$adj.matrix.fit = adj.matrix;
    probabilities.observed = list(marginal.probs=marginal.probs,joint.probs=joint.probs,conditional.probs=conditional.probs);
    probabilities.fit = list(estimated.marginal.probs=estimated.probabilities.fit$marginal.probs,estimated.joint.probs=estimated.probabilities.fit$joint.probs,estimated.conditional.probs=estimated.probabilities.fit$conditional.probs);
    probabilities = list(probabilities.observed=probabilities.observed,probabilities.fit=probabilities.fit);
    	
    # save the results for the model
    model[["caprese"]] = list(probabilities=probabilities,parents.pos=parents.pos,error.rates=error.rates,adj.matrix=adj.matrix.fit);
    
    # set the execution parameters
    parameters = list(algorithm="CAPRESE",lambda=lambda,do.estimation=do.estimation,silent=silent);
    
    #return the results
    topology = list(dataset=dataset,confidence=confidence,model=model,parameters=parameters,execution.time=(proc.time()-ptm));
    return(topology);
    
}

#' select at the most one parent for each node based on the probability raising criteria
#' @title get.tree.parents
#' @param adj.matrix adjacency matrix of the valid edges
#' @param marginal.probs observed marginal probabilities
#' @param joint.probs observed joint probabilities
#' @param lambda shrinkage parameter (value between 0 and 1)
#' @return best.parents list of the best parents
get.tree.parents = function( adj.matrix, marginal.probs, joint.probs, lambda ) {
	
    # compute the scores for each edge
    scores = get.tree.scores(adj.matrix,marginal.probs,joint.probs,lambda);
    pr.score = scores$pr.score;
    # set to -1 the scores where there is no causation according to Suppes' condition
    # [i,j] means i is causing j
    for (i in 1:ncol(pr.score)) {
        for (j in i:ncol(pr.score)) {
            # the diagonal has not to be considered (no self-cause)
            if(i==j) {
                pr.score[i,j] = -1;
            }
            # otherwise, apply Suppes's criteria for prima facie cause
            else {
                # if both the scores are not greater then 0, they are not valid
                # in this case the events are causally irrelevant, i.e., independent
                if(pr.score[i,j]<=0 && pr.score[j,i]<=0) {
                    pr.score[i,j] = -1;
                    pr.score[j,i] = -1;
                }
                # if at least one score is greater then 0, I keep the greater one
                # in this way I give a (time) direction to the progression
                # furthermore, this constrain the topology to be acyclic by construction
                else {
                    if(pr.score[i,j]>pr.score[j,i]) {
                        pr.score[j,i] = -1;
                    }
                    else {
                        pr.score[i,j] = -1;
                    }
                }
            }
        }
    }
    
    # chose at the most one parent per node
    # here I suppose that each node has a parent
    # spurious causes are considered (and removed) later
    best.parents = array(-1, dim=c(ncol(pr.score),1));
    for (i in 1:ncol(pr.score)) {
        # -1 means that the best parent is the Root
        curr.best = -1;
        # find the best parent for the current node
        best = which.max(pr.score[,i]);
        if(pr.score[best,i]>0) {
            curr.best = best;
        }
        # set the best parent for the current node
        best.parents[i,1] = curr.best;
    }
    
    # check for spurious causes by the independent progression filter and complete the parents list
    parents = verify.parents(best.parents,marginal.probs,joint.probs);
	
	# save the results    
    best.parents = list(parents=parents,marginal.probs=marginal.probs,joint.probs=joint.probs,pr.score=scores$pr.score);
    return(best.parents);
    
}

#' compute the probability raising based scores
#' @title get.tree.scores
#' @param adj.matrix adjacency matrix of the valid edges
#' @param marginal.probs observed marginal probabilities
#' @param joint.probs observed joint probabilities
#' @param lambda shrinkage parameter (value between 0 and 1)
#' @return scores: probability raising based scores
get.tree.scores = function( adj.matrix, marginal.probs, joint.probs, lambda ) {
	
    #structure where to save the probability raising scores
    pr.score = array(-1, dim=c(nrow(marginal.probs),nrow(marginal.probs)));
    
    # compute the probability raising based scores
    for (i in 1:ncol(pr.score)) {
        for (j in 1:ncol(pr.score)) {
        		
        		# if the edge is valid
        		if(adj.matrix[i,j]==1) {
        		
        			# alpha is the probability raising model of causation (raw model estimate)
        			alpha = ((joint.probs[i,j]/marginal.probs[i])-((marginal.probs[j]-joint.probs[i,j])/(1-marginal.probs[i])))/((joint.probs[i,j]/marginal.probs[i])+((marginal.probs[j]-joint.probs[i,j])/(1-marginal.probs[i])));
        			
        			# beta is the correction factor (based on time distance in terms of statistical dependence)
        			beta = (joint.probs[i,j]-marginal.probs[i]*marginal.probs[j])/(joint.probs[i,j]+marginal.probs[i]*marginal.probs[j]);
        			
        			# the overall estimator is a shrinkage-like combination of alpha and beta
        			# the scores are saved in the convention used for an ajacency matrix, i.e. [i,j] means causal edge i-->j
        			pr.score[i,j] = (1-lambda)*alpha + lambda*beta;
        			
        		}
            
        }
    }
    scores = list(marginal.probs=marginal.probs,joint.probs=joint.probs,pr.score=pr.score);
    return(scores);
    
}

#' verify the independent progression filter
#' @title verify.parents
#' @param best.parents best edges to be verified
#' @param marginal.probs observed marginal probabilities
#' @param joint.probs observed joint probabilities
#' @return best.parents: list of the best valid parents
verify.parents = function( best.parents, marginal.probs, joint.probs ) {
	
    # verify the condition for the best parent of each node
    for (i in 1:length(best.parents)) {
        # if there is a connection, i.e. the node is not already attached to the Root
        if(best.parents[i]!=-1) {
            # score for the root as the parent of this node
            w.root.node = 1/(1+marginal.probs[i]);
            # compute the scores for the edges to all the other upstream nodes
            attach.to.root = 1;
            for (j in 1:length(marginal.probs)) {
                # if the connection is valid and the parent node has greater probability
                # i.e. it is before the child in temporal order
                if(i!=j && marginal.probs[j]>marginal.probs[i]) {
                    w.parent.node = (marginal.probs[j]/(marginal.probs[i]+marginal.probs[j]))*(joint.probs[i,j]/(marginal.probs[i]*marginal.probs[j]));
                    # the parent is valid if this condition is valid at least one time (i.e. for at least one of the upstream nodes)
                    # meaning that if we find out that a connection is not spurious for any node, the best parent is not spurious as well
                    if(w.root.node<=w.parent.node) {
                        attach.to.root = 0;
                        break;
                    }
                }
            }
            # connect the node to the Root if the flag is true
            if(attach.to.root==1) {
                best.parents[i] = -1;
            }
        }
    }
    return(best.parents);
    
}
