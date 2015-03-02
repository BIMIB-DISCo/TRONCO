#### get.dag.scores.R
####
#### TRONCO: a tool for TRanslational ONCOlogy
####
#### See the files COPYING and LICENSE for copyright and licensing
#### information.


#compute the observed probabilities and the prima facie scores on the dataset
#INPUT:
#dataset: a valid dataset
#adj.matrix: adjacency matrix of the initially valid edges
#RETURN:
#scores: observed probabilities and prima facie scores
"get.dag.scores" <-
function(dataset, adj.matrix) {
	
    #structure to save the prima facie scores
    prima.facie.model <- array(-1, dim=c(ncol(dataset), ncol(dataset)));
    prima.facie.null <- array(-1, dim=c(ncol(dataset), ncol(dataset)));
    #structure to save the observed and observed-joint probabilities
    pair.count <- array(0, dim=c(ncol(dataset), ncol(dataset)));
    
	#compute the observed probabilities on the dataset
	for(i in 1:ncol(dataset)) {
	    for(j in 1:ncol(dataset)) {
	        val1 = dataset[ ,i];
	        val2 = dataset[ ,j];
	        pair.count[i,j] = (t(val1) %*% val2);
	    }
	}
    #marginal.probs is an array with the marginal probabilities
	marginal.probs <- array(as.matrix(diag(pair.count)/nrow(dataset)),dim=c(ncol(dataset),1));
    #joint.probs is an array with the joint observed probabilities
	joint.probs <- as.matrix(pair.count/nrow(dataset));
	
	#compute the prima facie scores based on the probability raising model
	for(i in 1:nrow(prima.facie.model)) {
	    for(j in 1:ncol(prima.facie.model)) {
            #the scores are saved in the convention of the adjacency matrix, i.e., [i,j] means i is causing j
            #the diagonal (self cause) and the other invalid edges have not to be considered
            if(adj.matrix[i,j]!=0) {
                #check if the connections from j to i and from i to j can be evaluated on this dataset
                if(marginal.probs[i]>0 && marginal.probs[i]<1 && marginal.probs[j]>0 && marginal.probs[j]<1) {
                    #check if the two events i and j are distinguishable
                    if((joint.probs[i,j]/marginal.probs[j])<1 || (joint.probs[i,j]/marginal.probs[i])<1) {
                        #prima facie scores of i --> j
                        prima.facie.model[i,j] = joint.probs[j,i]/marginal.probs[i];
                        prima.facie.null[i,j] = (marginal.probs[j]-joint.probs[j,i])/(1-marginal.probs[i]);
                    }
                }
            }
	    }
	}
	
	#save the results and return them
	scores <- list(marginal.probs=marginal.probs,joint.probs=joint.probs,prima.facie.model=prima.facie.model,prima.facie.null=prima.facie.null);
	return(scores);
	
}

#### end of file -- get.dag.scores.R
