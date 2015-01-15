#### get.tree.parents.R
####
#### TRONCO: a tool for TRanslational ONCOlogy
####
#### See the files COPYING and LICENSE for copyright and licensing
#### information.


#select at the most one parent for each node based on the probability raising criteria
#INPUT:
#marginal.probs: observed marginal probabilities
#joint.probs: observed joint probabilities
#lambda: shrinkage parameter (value between 0 and 1)
#RETURN:
#best.parents: list of the best parents
"get.tree.parents" <-
function(marginal.probs,joint.probs,lambda) {
	#compute the scores for each edge
	scores = get.tree.scores(marginal.probs,joint.probs,lambda);
	pr.score = scores$pr.score;
	#set to -1 the scores where there is no causation according to Suppes' condition
    #[i,j] means i is causing j
    for (i in 1:ncol(pr.score)) {
        for (j in i:ncol(pr.score)) {
            #the diagonal has not to be considered (no self-cause)
			if(i==j) {
				pr.score[i,j] = -1;
			}
            #otherwise, apply Suppes's criteria for prima facie cause
			else {
				#if both the scores are not greater then 0, they are not valid
                #in this case the events are causally irrelevant, i.e., independent
				if(pr.score[i,j]<=0 && pr.score[j,i]<=0) {
                    pr.score[i,j] = -1;
                    pr.score[j,i] = -1;
				}
				#if at least one score is greater then 0, I keep the greater one
                #in this way I give a (time) direction to the progression
                #furthermore, this constrain the topology to be acyclic by construction
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
	#chose at the most one parent per node
    #here I suppose that each node has a parent
    #spurious causes are considered (and removed) later
	best.parents = array(-1, dim=c(ncol(pr.score),1));
	for (i in 1:ncol(pr.score)) {
        #-1 means that the best parent is the Root
        curr.best = -1;
        #find the best parent for the current node
        best = which.max(pr.score[,i]);
        if(pr.score[best,i]>0) {
            curr.best = best;
        }
        #set the best parent for the current node
        best.parents[i,1] = curr.best;
	}
    #check for spurious causes by the independent progression filter and complete the parents list
    parents = verify.parents(best.parents,marginal.probs,joint.probs);
    best.parents = list(parents=parents,marginal.probs=marginal.probs,joint.probs=joint.probs,pr.score=scores$pr.score);
    return(best.parents);
}

#### end of file -- get.tree.parents.R
