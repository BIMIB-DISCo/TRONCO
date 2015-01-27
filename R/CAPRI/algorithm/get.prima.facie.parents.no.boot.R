#### get.prima.facie.parents.no.boot.R
####
#### TRONCO: a tool for TRanslational ONCOlogy
####
#### See the files COPYING and LICENSE for copyright and licensing
#### information.


#select the set of the prima facie parents (without bootstrap) for each node based on Suppes' definition of causation
#INPUT:
#dataset: a valid dataset
#adj.matrix: adjacency matrix of the initially valid edges
#RETURN:
#prima.facie.parents: list of the set (if any) of prima facie parents for each node
"get.prima.facie.parents.no.boot" <-
function(dataset, adj.matrix) {
	#compute the scores from the dataset
	scores = get.dag.scores(dataset,adj.matrix);
    #remove all the edges not representing a prima facie causes
    prima.facie.topology = get.prima.facie.causes.no.boot(adj.matrix,scores$marginal.probs,scores$prima.facie.model,scores$prima.facie.null);
    #save the results in a list and return it
    prima.facie.parents <- list(marginal.probs=scores$marginal.probs,joint.probs=scores$joint.probs,adj.matrix=prima.facie.topology,pf.confidence=NA);
    return(prima.facie.parents);
}

#### end of file -- get.prima.facie.parents.no.boot.R
