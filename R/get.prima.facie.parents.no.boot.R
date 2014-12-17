#### get.prima.facie.parents.no.boot.R
####
#### TRONCO: a tool for TRanslational ONCOlogy
####
#### See the files COPYING and LICENSE for copyright and licensing
#### information.


#select the set of the prima facie parents for each node based on Suppes' definition of causation without bootstrap
#INPUT:
#dataset: a valid dataset
#RETURN:
#prima.facie.parents: list of the set of prima facie parents for each node (if any)
"get.prima.facie.parents.no.boot" <-
function(dataset) {
	#compute the scores from the dataset
	scores = get.dag.scores(dataset);
    #remove all the edges not representing a prima facie causes
    prima.facie.topology = get.prima.facie.causes.no.boot(scores$marginal.probs,scores$prima.facie.model,scores$prima.facie.null);
    #save the results in a list and return it
    prima.facie.parents <- list(marginal.probs=scores$marginal.probs,joint.probs=scores$joint.probs,adj.matrix=prima.facie.topology,pf.confidence=NA);
    return(prima.facie.parents);
}

#### end of file -- get.prima.facie.parents.no.boot.R
