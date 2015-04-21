#### estimate.tree.error.rates.R
####
#### TRONCO: a tool for TRanslational ONCOlogy
####
#### See the files COPYING and LICENSE for copyright and licensing
#### information.


#estimate the error rates by "L-BFGS-B" optimization in terms of L2-error
#INPUT:
#marginal.probs: marginal probabilities
#joint.probs: joint probabilities
#parents.pos: which event is the parent? 0 if none, a number otherwise
#RETURN:
#estimated.error.rates: estimated probabilities, false positive and false negative error rates
"estimate.tree.error.rates" <-
function(marginal.probs,joint.probs,parents.pos) {
    #function to be optimized by "L-BFGS-B" optimization in terms of L2-error
	f.estimation <- function(errors) {
        #set the current error rates with the starting point of the optimization being (e_pos,e_neg) = 0
        error.rates = list(error.fp=errors[1],error.fn=errors[2]);
        #estimate the observed probabilities given the error rates
        estimated.probs = estimate.tree.probs(marginal.probs,joint.probs,parents.pos,error.rates);
        #evaluate the goodness of the estimatione by L2-error on the estimated marginal and joint probabilities
        error.estimation = sum((marginal.probs-estimated.probs$marginal.probs)^2)+sum((joint.probs-estimated.probs$joint.probs)^2);
        return(error.estimation);
	}
    #the estimation is performed as in Byrd et al (1995)
    #this method allows for box constraints, i.e., each variable can be given a lower and/or upper bound
	estimated.error.rates = optim(c(0.00,0.00),f.estimation,method="L-BFGS-B",lower=c(0.00,0.00),upper=c(0.49,0.49))$par;
    #structure to save the results
    estimated.error.rates = list(error.fp=estimated.error.rates[1],error.fn=estimated.error.rates[2]);
    return(estimated.error.rates);
}

#### end of file -- estimate.tree.error.rates.R


#### estimate.tree.joint.probs.R
####
#### TRONCO: a tool for TRanslational ONCOlogy
####
#### See the files COPYING and LICENSE for copyright and licensing
#### information.


#estimate the theoretical joint probability of two given nodes given the reconstructed topology
#INPUT:
#first.node: first node
#second.node: second node
#parents.pos: which event is the parent? -1 if none, a number otherwise
#marginal.probs: marginal probabilities
#conditional.probs: conditional probabilities
#RETURN:
#estimated.tree.joint.probs: estimated theoretical joint probability
"estimate.tree.joint.probs" <-
function(first.node,second.node,parents.pos,marginal.probs,conditional.probs) {
    #if the two nodes are roots
    if(parents.pos[first.node]==-1 && parents.pos[second.node]==-1) {
        estimated.tree.joint.probs = marginal.probs[first.node,1]*marginal.probs[second.node,1];
    }
    #if the two nodes are the same node
    else if(first.node==second.node) {
        estimated.tree.joint.probs = marginal.probs[first.node,1];
    }
    #otherwise
    else {
        #go through the parents starting from the two nodes to find if they are directly connected
        #are the two nodes in the same path?
        is.path = 0;
        #check if first.node is an ancestor of second.node
        curr.first = first.node;
        curr.second = second.node;
        while(parents.pos[curr.second]!=-1) {
            if(curr.first==curr.second) {
                is.path = 1;
                is.child = second.node;
                break;
            }
            curr.second = parents.pos[curr.second];
        }
        if(is.path==0) {
            #check if second.node is an ancestor of first.node
            curr.first = first.node;
            curr.second = second.node;
            while(parents.pos[curr.first]!=-1) {
                if(curr.first==curr.second) {
                    is.path = 1;
                    is.child = first.node;
                    break;
                }
                curr.first = parents.pos[curr.first];
            }
        }
        #check if the two nodes are at least connected
        #are the two nodes connected?
        is.connected = 0;
        if(is.path==0) {
            curr.first = first.node;
            curr.second = second.node;
            while(parents.pos[curr.first]!=-1) {
                while(parents.pos[curr.second]!=-1) {
                    if(curr.first==curr.second) {
                        is.connected = 1;
                        is.ancestor = curr.first;
                        break;
                    }
                    else {
                        curr.second = parents.pos[curr.second];
                    }
                }
                if(is.connected==0) {
                    curr.first = parents.pos[curr.first];
                    curr.second = second.node;
                }
                else {
                    break;
                }
            }
        }
        #now I can set the joint probabilities
        #in this case the two nodes are directly connected
        #P(child,parent)_estimate = P(child);
        if(is.path==1) {
            estimated.tree.joint.probs = marginal.probs[is.child,1];
        }
        #in this case the two nodes are indirectly connected
        #P(i,j)_estimate = P(ancestor)_estimate * P_PATH(ancestor->first.node)_estimate * P_PATH(ancestor->second.node)_estimate
        else if(is.connected==1) {
            #P(ancestor)_estimate
            estimated.tree.joint.probs = marginal.probs[is.ancestor,1];
            #P_PATH(ancestor->first.node)_estimate
            first.path = 1;
            curr.first = first.node;
            while(parents.pos[curr.first]!=is.ancestor) {
                first.path = first.path * conditional.probs[curr.first,1];
                curr.first = parents.pos[curr.first];
            }
            #P_PATH(ancestor->second.node)_estimate
            second.path = 1;
            curr.second = second.node;
            while(parents.pos[curr.second]!=is.ancestor) {
                second.path = second.path * conditional.probs[curr.second,1];
                curr.second = parents.pos[curr.second];
            }
            estimated.tree.joint.probs = estimated.tree.joint.probs * first.path * second.path;
        }
        #in this case the two nodes are not connected
        #P(i,j)_estimate = P(i)_estimate * P(j)_estimate
        else {
            estimated.tree.joint.probs = marginal.probs[first.node,1]*marginal.probs[second.node,1];
        }
    }
    return(estimated.tree.joint.probs);
}

#### end of file -- estimate.tree.joint.probs.R


#### estimate.tree.probs.R
####
#### TRONCO: a tool for TRanslational ONCOlogy
####
#### See the files COPYING and LICENSE for copyright and licensing
#### information.


#estimate the marginal, joint and conditional probabilities given the reconstructed topology and the error rates
#INPUT:
#marginal.probs: observed marginal probabilities
#joint.probs: observed joint probabilities
#parents.pos: position of the parents in the list of nodes
#error.rates: rates for the false positive and the false negative errors
#RETURN:
#estimated.probs: estimated marginal, joint and conditional probabilities
"estimate.tree.probs" <-
function(marginal.probs,joint.probs,parents.pos,error.rates) {
    #structure where to save the probabilities to be estimated
    estimated.marginal.probs = array(-1, dim=c(nrow(marginal.probs),1));
    estimated.joint.probs = array(-1, dim=c(nrow(marginal.probs),nrow(marginal.probs)));
    estimated.conditional.probs = array(-1, dim=c(nrow(marginal.probs),1));
    #estimate the theoretical conditional probabilities given the error rates
    #this estimation is performed by applying the error rates to the marginal and joint probabilities
    theoretical.conditional.probs = array(-1, dim=c(nrow(marginal.probs),1));
    for (i in 1:nrow(theoretical.conditional.probs)) {
        #if the node has a parent, use the error rates to compute the conditional probability
        #if the node has no parent, its conditional probability is not considered
        if(parents.pos[i,1]!=-1) {
            #P(i|j)_theoretical = ((P(i,j)_obs-e_p*(P(j)_obs+P(i)_obs)+e_p^2)/(1-e_n-e_p)^2)/((P(j)_obs-e_p)/(1-e_n-e_p))
            theoretical.conditional.probs[i,1] = (joint.probs[i,parents.pos[i,1]]-error.rates$error.fp*(marginal.probs[parents.pos[i,1],1]+marginal.probs[i,1])+error.rates$error.fp^2)/((marginal.probs[parents.pos[i,1],1]-error.rates$error.fp)*(1-error.rates$error.fn-error.rates$error.fp));
            if(theoretical.conditional.probs[i,1]<0 || theoretical.conditional.probs[i,1]>1) {
                #invalid theoretical conditional probability
                if(theoretical.conditional.probs[i,1]<0) {
                    theoretical.conditional.probs[i,1] = 0;
                }
                else {
                    theoretical.conditional.probs[i,1] = 1;
                }
            }
        }
    }
    #estimate the marginal observed probabilities
    #this estimation is performed by applying the topological constraints on the probabilities and then the error rates
    #I do not have any constraint on the nodes without a parent
    child.list <- which(parents.pos==-1);
    estimated.marginal.probs[child.list,1] = marginal.probs[child.list,1];
    estimated.marginal.probs.with.error = array(-1, dim=c(nrow(marginal.probs),1));
    estimated.marginal.probs.with.error[child.list,1] = estimated.marginal.probs[child.list,1];
    visited = length(child.list);
    #I do not have any constraint for the joint probabilities on the pair of nodes which are the roots of the tree/forest
    estimated.joint = array(0, dim=c(nrow(marginal.probs),nrow(marginal.probs)));
    for (i in child.list) {
        for (j in child.list) {
            if(i!=j) {
                estimated.joint.probs[i,j] = joint.probs[i,j];
                estimated.joint[i,j] = -1;
            }
        }
    }
    #visit the nodes with a parent in topological order
    while (visited < nrow(estimated.marginal.probs)) {
        #set the new child list
        new.child = vector();
        #go through the current parents
        for (node in child.list) {
            #set the new children
            curr.child <- which(parents.pos==node);
            #go through the current children
            for (child in curr.child) {
                #set the marginal probability for this node
                #P(child)_estimate = P(parent)_estimate * P(child|parent)_theoretical
                estimated.marginal.probs[child,1] = estimated.marginal.probs[parents.pos[child,1],1]*theoretical.conditional.probs[child,1];
                visited = visited + 1;
                #P(child,parent)_estimare = P(child)_estimate;
                estimated.joint.probs[child,parents.pos[child,1]] = estimated.marginal.probs[child,1];
                estimated.joint[child,parents.pos[child,1]] = 1;
                estimated.joint.probs[parents.pos[child,1],child] = estimated.marginal.probs[child,1];
                estimated.joint[parents.pos[child,1],child] = 1;
                #apply the error rates to the marginal probabilities
                #P(i)_obs_estimate = P(i)_estimate*(1-e_n) + P(not i)_estimate*e_p
                estimated.marginal.probs.with.error[child,1] = error.rates$error.fp+(1-error.rates$error.fn-error.rates$error.fp)*estimated.marginal.probs[child,1];
                if(estimated.marginal.probs.with.error[child,1]<0 || estimated.marginal.probs.with.error[child,1]>1) {
                    #invalid estimated observed probability
                    if(estimated.marginal.probs.with.error[child,1]<0) {
                        estimated.marginal.probs.with.error[child,1] = 0;
                    }
                    else {
                        estimated.marginal.probs.with.error[child,1] = 1;
                    }
                }
            }
            new.child <- c(new.child,curr.child);
        }
        #set the next child list
        child.list = new.child;
    }
    diag(estimated.joint.probs) = estimated.marginal.probs;
    diag(estimated.joint) = -1;
    #given the estimated observed probabilities, I can now also estimate the joint probabilities by applying the topological constraints and then the error rates
    for (i in 1:nrow(estimated.joint.probs)) {
        for (j in i:nrow(estimated.joint.probs)) {
            #if I still need to estimate this joint probability
            if(estimated.joint[i,j]==0) {
                estimated.joint.probs[i,j] = estimate.tree.joint.probs(i,j,parents.pos,estimated.marginal.probs,theoretical.conditional.probs);
                estimated.joint[i,j] = 1;
            }
            #now I can apply the error rates to estimate the observed joint probabilities
            if(estimated.joint[i,j]==1) {
                #P(i,j)_obs_estimate = P(i,j)_estimate*(1-e_n)^2+P(not i,j)_estimate*e_p*(1-e_n)+P(i,not j)_estimate*(1-e_n)*e_p+P(not i,not j)_estimate*e_p^2;
                estimated.joint.probs[i,j] = estimated.joint.probs[i,j]*((1-error.rates$error.fn-error.rates$error.fp)^2)+error.rates$error.fp*(estimated.marginal.probs[i,1]+estimated.marginal.probs[j,1])-error.rates$error.fp^2;
                #invalid estimated joint probability
                if(estimated.joint.probs[i,j]<0 || estimated.joint.probs[i,j]>min(estimated.marginal.probs.with.error[i,1],estimated.marginal.probs.with.error[j,1])) {
                    if(estimated.joint.probs[i,j]<0) {
                        estimated.joint.probs[i,j] = 0;
                    }
                    else {
                        estimated.joint.probs[i,j] = min(estimated.marginal.probs.with.error[i,1],estimated.marginal.probs.with.error[j,1]);
                    }
                }
                estimated.joint.probs[j,i] = estimated.joint.probs[i,j];
            }
        }
    }
    #save the estimated probabilities
    estimated.marginal.probs = estimated.marginal.probs.with.error;
    #given the estimated observed and joint probabilities, I can finally compute the conditional probabilities
    #P(child|parent)_obs_estimate = P(parent,child)_obs_estimate/P(parent)_obs_estimate
    for (i in 1:nrow(estimated.conditional.probs)) {
        if(parents.pos[i,1]!=-1) {
            if(estimated.marginal.probs[parents.pos[i,1],1]>0) {
                estimated.conditional.probs[i,1] = estimated.joint.probs[parents.pos[i,1],i]/estimated.marginal.probs[parents.pos[i,1],1];
            }
            else {
                estimated.conditional.probs[i,1] = 0;
            }
        }
        #if the node has no parent, its conditional probability is set to 1
        else {
            estimated.conditional.probs[i,1] = 1;
        }
    }
    #structure to save the results
    estimated.probs = list(marginal.probs=estimated.marginal.probs,joint.probs=estimated.joint.probs,conditional.probs=estimated.conditional.probs);
    return(estimated.probs);
}

#### end of file -- estimate.tree.probs.R


