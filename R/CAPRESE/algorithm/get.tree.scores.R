#### get.tree.scores.R
####
#### TRONCO: a tool for TRanslational ONCOlogy
####
#### See the files COPYING and LICENSE for copyright and licensing
#### information.


#compute the probability raising based scores
#INPUT:
#marginal.probs: observed marginal probabilities
#joint.probs: observed joint probabilities
#lambda: shrinkage parameter (value between 0 and 1)
#RETURN:
#scores: probability raising based scores
"get.tree.scores" <-
function(marginal.probs,joint.probs,lambda) {
    #structure where to save the probability raising scores
	pr.score = array(-1, dim=c(nrow(marginal.probs),nrow(marginal.probs)));
	#compute the probability raising based scores
	for (i in 1:ncol(pr.score)) {
        for (j in 1:ncol(pr.score)) {
            #alpha is the probability raising model of causation (raw model estimate)
            alpha = ((joint.probs[i,j]/marginal.probs[i])-((marginal.probs[j]-joint.probs[i,j])/(1-marginal.probs[i])))/((joint.probs[i,j]/marginal.probs[i])+((marginal.probs[j]-joint.probs[i,j])/(1-marginal.probs[i])));
            #beta is the correction factor (based on time distance in terms of statistical dependence)
            beta = (joint.probs[i,j]-marginal.probs[i]*marginal.probs[j])/(joint.probs[i,j]+marginal.probs[i]*marginal.probs[j]);
            #the overall estimator is a shrinkage-like combination of alpha and beta
            #the scores are saved in the convention used for an ajacency matrix, i.e. [i,j] means causal edge i-->j
            pr.score[i,j] = (1-lambda)*alpha + lambda*beta;
        }
	}
	scores = list(marginal.probs=marginal.probs,joint.probs=joint.probs,pr.score=pr.score);
	return(scores);
}

#### end of file -- get.tree.scores.R
