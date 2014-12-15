#### estimate.dag.joint.probs.R
####
#### TRONCO: a tool for TRanslational ONCOlogy
####
#### See the files COPYING and LICENSE for copyright and licensing
#### information.


#estimate the theoretical joint probability of two given nodes given the reconstructed topology
#INPUT:
#first.node: first node
#second.node: second node
#parents.pos: which event is the parent? -1 if none, a list otherwise
#marginal.probs: marginal probabilities
#conditional.probs: conditional probabilities
#RETURN:
#estimated.dag.joint.probs: estimated theoretical joint probability
"estimate.dag.joint.probs" <-
function(first.node,second.node,parents.pos,marginal.probs,conditional.probs) {
    #if the two nodes are roots
    if((length(parents.pos[[first.node]])==1 && parents.pos[[first.node]]==-1) && (length(parents.pos[[second.node]])==1 && parents.pos[[second.node]]==-1)) {
        estimated.dag.joint.probs = marginal.probs[first.node,1]*marginal.probs[second.node,1];
    }
    #if the two nodes are the same node
    else if(first.node==second.node) {
        estimated.dag.joint.probs = marginal.probs[first.node,1];
    }
    #otherwise
    else {
        #go through the parents starting from the two nodes to find if they are directly connected
        #are the two nodes in the same path?
        is.path = 0;
        #check if first.node is an ancestor of second.node
        curr.first = first.node;
        curr.second = second.node;
        while(length(unlist(parents.pos[curr.second]))>0 && (length(unlist(parents.pos[curr.second]))>1 || unlist(parents.pos[curr.second])!=-1)) {
			curr.result = which(curr.first%in%curr.second);
            if(length(curr.result)>0) {
                is.path = 1;
                is.child = second.node;
                break;
            }
            curr.second = unique(unlist(parents.pos[curr.second]));
            curr.second = curr.second[which(curr.second!=-1)];
        }
        if(is.path==0) {
            #check if second.node is an ancestor of first.node
            curr.first = first.node;
            curr.second = second.node;
            while(length(unlist(parents.pos[curr.first]))>0 && (length(unlist(parents.pos[curr.first]))>1 || unlist(parents.pos[curr.first])!=-1)) {
				curr.result = which(curr.first%in%curr.second);
                if(length(curr.result)>0) {
                    is.path = 1;
                    is.child = first.node;
                    break;
                }
				curr.first = unique(unlist(parents.pos[curr.first]));
	            curr.first = curr.first[which(curr.first!=-1)];
            }
        }
        #check if the two nodes are at least connected
        #are the two nodes connected?
        is.connected = 0;
        if(is.path==0) {
            curr.first = first.node;
            curr.second = second.node;
            while(length(unlist(parents.pos[curr.first]))>0 && (length(unlist(parents.pos[curr.first]))>1 || unlist(parents.pos[curr.first])!=-1)) {
                while(length(unlist(parents.pos[curr.second]))>0 && (length(unlist(parents.pos[curr.second]))>1 || unlist(parents.pos[curr.second])!=-1)) {
					curr.result = which(curr.first%in%curr.second);
                    if(length(curr.result)>0) {
                        is.connected = 1;
                        is.ancestor = curr.result;
                        break;
                    }
                    else {
						curr.second = unique(unlist(parents.pos[curr.second]));
						curr.second = curr.second[which(curr.second!=-1)];
                    }
                }
                if(is.connected==0) {
            		curr.first = unique(unlist(parents.pos[curr.first]));
            		curr.first = curr.first[which(curr.first!=-1)];
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
            estimated.dag.joint.probs = marginal.probs[is.child,1];
        }
        #in this case the two nodes are indirectly connected
        #P(i,j)_estimate = P(ancestor)_estimate * P_PATH(ancestor->first.node)_estimate * P_PATH(ancestor->second.node)_estimate
        else if(is.connected==1) {
			#as heuristic, if I have multiple common ancestors, I choose the later one to estimate the joint probability
			is.ancestor = is.ancestor[which.min(marginal.probs[is.ancestor,1])];
            #P(ancestor)_estimate
            estimated.dag.joint.probs = marginal.probs[is.ancestor,1];
            #P_PATH(ancestor->first.node)_estimate
            all.first.paths = enumerate.all.paths(is.ancestor,first.node,parents.pos);
            first.path = rep(1,length(all.first.paths));
            for (i in 1:length(all.first.paths)) {
				curr.new.path = all.first.paths[[i]];
				if(length(curr.new.path)>0) {
					for (j in 2:length(curr.new.path)) {
						curr.parents = parents.pos[[curr.new.path[j-1]]];
						curr.conditional.probs = conditional.probs[[curr.new.path[j],1]];
						first.path[i] = first.path[i] * curr.conditional.probs;
					}
				}
            }
            #as heuristic, if I have multiple paths, I average the probability of each path
            first.path = mean(first.path);
            #P_PATH(ancestor->second.node)_estimate
            all.second.paths = enumerate.all.paths(is.ancestor,second.node,parents.pos);
            second.path = rep(1,length(all.second.paths));
            for (i in 1:length(all.second.paths)) {
				curr.new.path = all.second.paths[[i]];
				if(length(curr.new.path)>0) {
					for (j in 2:length(curr.new.path)) {
						curr.parents = parents.pos[[curr.new.path[j-1]]];
						curr.conditional.probs = conditional.probs[[curr.new.path[j],1]];
						second.path[i] = second.path[i] * curr.conditional.probs;
					}
				}
            }
            #as heuristic, if I have multiple paths, I average the probability of each path
            second.path = mean(second.path);
            estimated.dag.joint.probs = estimated.dag.joint.probs * first.path * second.path;
        }
        #in this case the two nodes are not connected
        #P(i,j)_estimate = P(i)_estimate * P(j)_estimate
        else {
            estimated.dag.joint.probs = marginal.probs[first.node,1]*marginal.probs[second.node,1];
        }
    }
    return(estimated.dag.joint.probs);
}

#### end of file -- estimate.dag.joint.probs.R
