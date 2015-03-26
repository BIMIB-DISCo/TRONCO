#### get.prima.facie.parents.R
####
#### TRONCO: a tool for TRanslational ONCOlogy
####
#### See the files COPYING and LICENSE for copyright and licensing
#### information.


#select the set of the prima facie parents (with bootstrap) for each node based on Suppes' definition of causation
#INPUT:
#dataset: a valid dataset
#nboot: integer number (greater than 0) of bootstrap sampling to be performed
#pvalue: pvalue for the tests (value between 0 and 1)
#adj.matrix: adjacency matrix of the initially valid edges
#RETURN:
#prima.facie.parents: list of the set (if any) of prima facie parents for each node
"get.prima.facie.parents.do.boot" <-
function(dataset, hypotheses, nboot, pvalue, adj.matrix) {
	
	#compute a robust estimation of the scores using rejection sampling bootstrap
	scores = get.bootstapped.scores(dataset,nboot,adj.matrix);
	
    #remove all the edges not representing a prima facie causes
    prima.facie.topology = get.prima.facie.causes.do.boot(adj.matrix,hypotheses,scores$marginal.probs.distributions,scores$prima.facie.model.distributions,scores$prima.facie.null.distributions,pvalue);
    
    #compute the observed and joint probabilities as the mean of the bootstrapped values
    marginal.probs = array(-1,dim=c(ncol(dataset),1));
    joint.probs = array(-1,dim=c(ncol(dataset),ncol(dataset)));
    for(i in 1:ncol(dataset)) {
        marginal.probs[i,1] = mean(unlist(scores$marginal.probs.distributions[i,1]));
        for(j in i:ncol(dataset)) {
            joint.probs[i,j] = mean(unlist(scores$joint.probs.distributions[i,j]));
            if(i!=j) {
                joint.probs[j,i] = joint.probs[i,j];
            }
        }
    }
    
    #save the results and return them
    prima.facie.parents <- list(marginal.probs=marginal.probs,joint.probs=joint.probs,adj.matrix=prima.facie.topology$adj.matrix,pf.confidence=prima.facie.topology$edge.confidence.matrix);
    return(prima.facie.parents);
    
}


#select the set of the prima facie parents (without bootstrap) for each node based on Suppes' definition of causation
#INPUT:
#dataset: a valid dataset
#adj.matrix: adjacency matrix of the initially valid edges
#RETURN:
#prima.facie.parents: list of the set (if any) of prima facie parents for each node
"get.prima.facie.parents.no.boot" <-
function(dataset, hypotheses, adj.matrix) {
	
	#compute the scores from the dataset
	scores = get.dag.scores(dataset,adj.matrix);
	
    #remove all the edges not representing a prima facie causes
    prima.facie.topology = get.prima.facie.causes.no.boot(adj.matrix,hypotheses,scores$marginal.probs,scores$prima.facie.model,scores$prima.facie.null);
    
    #save the results return them
    prima.facie.parents <- list(marginal.probs=scores$marginal.probs,joint.probs=scores$joint.probs,adj.matrix=prima.facie.topology,pf.confidence=NA);
    return(prima.facie.parents);
    
}

#### end of file -- get.prima.facie.parents.R
