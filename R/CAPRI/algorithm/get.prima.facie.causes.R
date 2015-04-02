#### get.prima.facie.causes.R
####
#### TRONCO: a tool for TRanslational ONCOlogy
####
#### See the files COPYING and LICENSE for copyright and licensing
#### information.


# select the best set of prima facie causes per node
# INPUT:
# adj.matrix: adjacency matrix of the initially valid edges
# hypotheses: hypotheses to be considered
# marginal.probs.distributions: distributions of the bootstrapped marginal probabilities
# prima.facie.model.distributions: distributions of the prima facie model
# prima.facie.null.distributions: distributions of the prima facie null
# pvalue: minimum pvalue for the Mann-Whitney U tests to be significant
# dataset: a valid dataset
# marginal.probs: observed marginal probabilities
# joint.probs: observed joint probabilities
# RETURN:
# prima.facie.topology: list describing the topology of the prima facie causes
"get.prima.facie.causes.do.boot" <-
function( adj.matrix, hypotheses, marginal.probs.distributions, prima.facie.model.distributions, prima.facie.null.distributions, pvalue, dataset, marginal.probs, joint.probs, silent = FALSE ) {
	
    # structure to save the confidence of the edges
    edge.confidence.matrix <- array(list(), c(3,1));
    edge.confidence.matrix[[1,1]] = array(NA, c(ncol(prima.facie.model.distributions),ncol(prima.facie.model.distributions)));
    edge.confidence.matrix[[2,1]] = array(NA, c(ncol(prima.facie.model.distributions),ncol(prima.facie.model.distributions)));
    edge.confidence.matrix[[3,1]] = array(NA, c(ncol(prima.facie.model.distributions),ncol(prima.facie.model.distributions)));
    
    # verify Suppes' conditions for prima facie causes
    # i.e., i --> j implies P(i)>P(j) (temporal priority) and P(j|i)>P(j|not i) (probability raising)
    # verify the temporal priority condition
    if(!silent) cat(paste0('\tEvaluating \"temporal priority\" (Wilcoxon, p-value ', pvalue, ')\n'));
    temporal.priority = verify.temporal.priority.do.boot(marginal.probs.distributions,pvalue,adj.matrix,edge.confidence.matrix);
    
    # verify the probability raising condition
    if(!silent) cat(paste0('\tEvaluating \"probability raising\" (Wilcoxon, p-value ', pvalue, ')\n'));
    probability.raising = verify.probability.raising.do.boot(prima.facie.model.distributions,prima.facie.null.distributions,pvalue,temporal.priority$adj.matrix,temporal.priority$edge.confidence.matrix);
    
    # perform the hypergeometric test for each pair of events
    for(i in 1:ncol(adj.matrix)) {
		for(j in i:nrow(adj.matrix)) {
			
			# the diagonal (self cause) has not to be considered
			if(i!=j) {
				#compute the confidence by hypergeometric test for both j --> i and i --> j
				probability.raising$edge.confidence.matrix[[3,1]][i,j] = phyper(joint.probs[i,j]*nrow(dataset),marginal.probs[i]*nrow(dataset),nrow(dataset)-marginal.probs[i]*nrow(dataset),marginal.probs[j]*nrow(dataset),lower.tail=FALSE);
				probability.raising$edge.confidence.matrix[[3,1]][j,i] = probability.raising$edge.confidence.matrix[[3,1]][i,j];
			}
        
		}
    }
    
    # remove any cycle
    if(length(temporal.priority$not.ordered)>0 || !is.na(hypotheses[1])) {
    		if(!silent) cat('*** Loop detection found loops to break.\n')
    		weights.temporal.priority = probability.raising$edge.confidence.matrix[[1,1]]+probability.raising$edge.confidence.matrix[[2,1]];
        weights.matrix = probability.raising$edge.confidence.matrix[[2,1]]+probability.raising$edge.confidence.matrix[[3,1]];
        acyclic.topology = remove.cycles(probability.raising$adj.matrix,weights.temporal.priority,weights.matrix,temporal.priority$not.ordered,hypotheses);
        adj.matrix = acyclic.topology$adj.matrix;
    }
    else {
        adj.matrix = probability.raising$adj.matrix;
    }
    
    # save the results and return them
    prima.facie.topology <- list(adj.matrix=adj.matrix,edge.confidence.matrix=probability.raising$edge.confidence.matrix);
    return(prima.facie.topology);

}


# select the best set of prima facie causes per node without bootstrap
# INPUT:
# adj.matrix: adjacency matrix of the initially valid edges
# marginal.probs: marginal probabilities
# prima.facie.model: prima facie model
# prima.facie.null: prima facie null
# dataset: a valid dataset
# marginal.probs: observed marginal probabilities
# joint.probs: observed joint probabilities
# RETURN:
# prima.facie.topology: adjacency matrix of the prima facie causes
"get.prima.facie.causes.no.boot" <-
function( adj.matrix, hypotheses, marginal.probs, prima.facie.model, prima.facie.null, dataset, joint.probs, silent = FALSE ) {
	
    # structure to save the confidence of the edges
    edge.confidence.matrix <- array(list(), c(3,1));
    edge.confidence.matrix[[1,1]] = array(NA, c(ncol(adj.matrix),ncol(adj.matrix)));
    edge.confidence.matrix[[2,1]] = array(NA, c(ncol(adj.matrix),ncol(adj.matrix)));
    edge.confidence.matrix[[3,1]] = array(NA, c(ncol(adj.matrix),ncol(adj.matrix)));
    
    # verify Suppes' conditions for prima facie causes
    # i.e., i --> j implies P(i)>P(j) (temporal priority) and P(j|i)>P(j|not i) (probability raising)
    # verify the temporal priority condition
    if(!silent) cat(paste0('\tEvaluating \"temporal priority\".\n'));
    temporal.priority = verify.temporal.priority.no.boot(marginal.probs,adj.matrix,edge.confidence.matrix);
    
    # verify the probability raising condition
    if(!silent) cat(paste0('\tEvaluating \"probability raising\".\n'));
    probability.raising = verify.probability.raising.no.boot(prima.facie.model,prima.facie.null,temporal.priority$adj.matrix,temporal.priority$edge.confidence.matrix);
    
    # perform the hypergeometric test for each pair of events
    for(i in 1:ncol(adj.matrix)) {
		for(j in i:nrow(adj.matrix)) {
			
			# the diagonal (self cause) has not to be considered
			if(i!=j) {
				#compute the confidence by hypergeometric test for both j --> i and i --> j
				probability.raising$edge.confidence.matrix[[3,1]][i,j] = phyper(joint.probs[i,j]*nrow(dataset),marginal.probs[i]*nrow(dataset),nrow(dataset)-marginal.probs[i]*nrow(dataset),marginal.probs[j]*nrow(dataset),lower.tail=FALSE);
				probability.raising$edge.confidence.matrix[[3,1]][j,i] = probability.raising$edge.confidence.matrix[[3,1]][i,j];
			}
        
		}
    }
    
    # remove any cycle
    if(length(temporal.priority$not.ordered)>0 || !is.na(hypotheses[1])) {
    		if(!silent) cat('*** Loop detection found loops to break.\n')
    		weights.temporal.priority = probability.raising$edge.confidence.matrix[[2,1]];
        weights.matrix = probability.raising$edge.confidence.matrix[[3,1]];
        acyclic.topology = remove.cycles(probability.raising$adj.matrix,weights.temporal.priority,weights.matrix,temporal.priority$not.ordered,hypotheses);
        adj.matrix = acyclic.topology$adj.matrix;
    }
    else {
        adj.matrix = probability.raising$adj.matrix;
    }
    
    # save the results and return them
    prima.facie.topology <- list(adj.matrix=adj.matrix,edge.confidence.matrix=edge.confidence.matrix);
    return(prima.facie.topology);

}

#### end of file -- get.prima.facie.causes.R
