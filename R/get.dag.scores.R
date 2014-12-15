#### get.dag.scores.R
####
#### TRONCO: a tool for TRanslational ONCOlogy
####
#### See the files COPYING and LICENSE for copyright and licensing
#### information.


#compute the observed probabilities and the prima facie scores on the dataset
#INPUT:
#dataset: a valid dataset
#RETURN:
#scores: observed probabilities and prima facie scores
"get.dag.scores" <-
function(dataset) {
    #structure to save the positive prima facie scores
    prima.facie.model <- array(-1, dim=c(ncol(dataset), ncol(dataset)));
    prima.facie.null <- array(-1, dim=c(ncol(dataset), ncol(dataset)));
    #structure to compute the observed and observed-joint probabilities
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
	#compute the prima facie errors based on the probability raising model
	for(i in 1:nrow(prima.facie.model)) {
	    for(j in i:ncol(prima.facie.model)) {
            #the scores are saved in the convention of the adjacency matrix, i.e. [j,i] means j is causing i
            #the diagonal (self cause) has not to be considered
            if(i!=j) {
                #check if the connections from j to i and from i to j can be evaluated on this dataset
                if(marginal.probs[i]>0 && marginal.probs[i]<1 && marginal.probs[j]>0 && marginal.probs[j]<1) {
                    #check if the two events i and j are distinguishable
                    if((joint.probs[i,j]/marginal.probs[j])<1 || (joint.probs[i,j]/marginal.probs[i])<1) {
                        #prima facie scores of i --> j
                        prima.facie.model[i,j] = joint.probs[j,i]/marginal.probs[i];
                        prima.facie.null[i,j] = (marginal.probs[j]-joint.probs[j,i])/(1-marginal.probs[i]);
                        #prima facie scores of j --> i
                        prima.facie.model[j,i] = joint.probs[i,j]/marginal.probs[j];
                        prima.facie.null[j,i] = (marginal.probs[i]-joint.probs[i,j])/(1-marginal.probs[j]);
                    }
                }
            }
            else {
                prima.facie.model[i,j] = 1;
                prima.facie.null[i,j] = 0;
            }
	    }
	}
	scores <- list(marginal.probs=marginal.probs,joint.probs=joint.probs,prima.facie.model=prima.facie.model,prima.facie.null=prima.facie.null);
	return(scores);
}

#### end of file -- get.dag.scores.R
