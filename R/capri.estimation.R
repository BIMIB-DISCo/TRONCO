#### enumerate.all.paths.R
####
#### TRONCO: a tool for TRanslational ONCOlogy
####
#### See the files COPYING and LICENSE for copyright and licensing
#### information.


#enumerate all the paths between two nodes of a DAG
#INPUT:
#ancestor.node: first node of the path
#child.node: last node of the path
#parents.pos: topological connections
#RETURN:
#all.paths: vector of all the paths
"enumerate.all.paths" <-
function( ancestor.node, child.node, parents.pos ) {
	#set the initial parents set
	all.paths = parents.pos[[child.node]];
	if(length(all.paths)==1 && all.paths==-1) {
		is.done = TRUE;
	}
	else {
		is.done = FALSE;
	}
	#visit all the nodes in topological order
	while(is.done==FALSE) {
		curr.paths = vector();
		is.done = TRUE;
		for (i in 1:length(all.paths)) {
			curr.new.path = all.paths[[i]];
			curr.new.parents = parents.pos[[curr.new.path[1]]];
			if(length(curr.new.parents)>1 || curr.new.parents!=-1) {
				is.done = FALSE;
				for (j in 1:length(curr.new.parents)) {
					curr.paths[length(curr.paths)+1] = list(c(curr.new.parents[j],curr.new.path));
				}
			}
			else {
				curr.paths[length(curr.paths)+1] = list(curr.new.path);
			}
		}
		all.paths = curr.paths;
	}
	#remove all the paths that are not visiting ancestor.node
	curr.paths = vector();
	for (i in 1:length(all.paths)) {
		curr.result = which(all.paths[[i]]%in%ancestor.node);
		if(length(curr.result)>0) {
			curr.new.path = all.paths[[i]];
			curr.paths[length(curr.paths)+1] = list(c(curr.new.path[which(all.paths[[i]]%in%ancestor.node):length(curr.new.path)],child.node));
		}
	}
	all.paths = unique(curr.paths);
    return(all.paths);
}

#### end of file -- enumerate.all.paths.R


#### estimate.dag.error.rates.R
####
#### TRONCO: a tool for TRanslational ONCOlogy
####
#### See the files COPYING and LICENSE for copyright and licensing
#### information.


#estimate the error rates by "L-BFGS-B" optimization in terms of L2-error
#INPUT:
#dataset: a valid dataset
#marginal.probs: marginal probabilities
#joint.probs: joint probabilities
#parents.pos: which event is the parent? 0 if none, a number otherwise
#RETURN:
#estimated.error.rates: estimated probabilities, false positive and false negative error rates
"estimate.dag.error.rates" <-
function( dataset, marginal.probs, joint.probs, parents.pos ) {
    #function to be optimized by "L-BFGS-B" optimization in terms of L2-error
	f.estimation <- function(errors) {
        #set the current error rates with the starting point of the optimization being (e_pos,e_neg) = 0
        error.rates = list(error.fp=errors[1],error.fn=errors[2]);
        #estimate the observed probabilities given the error rates
        estimated.probs = estimate.dag.probs(dataset,marginal.probs,joint.probs,parents.pos,error.rates);
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

#### end of file -- estimate.dag.error.rates.R


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
function( first.node, second.node, parents.pos, marginal.probs, conditional.probs ) {
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


#### estimate.dag.probs.R
####
#### TRONCO: a tool for TRanslational ONCOlogy
####
#### See the files COPYING and LICENSE for copyright and licensing
#### information.


#estimate the marginal, joint and conditional probabilities given the reconstructed topology and the error rates
#INPUT:
#dataset: a valid dataset
#marginal.probs: observed marginal probabilities
#joint.probs: observed joint probabilities
#parents.pos: position of the parents in the list of nodes
#error.rates: rates for the false positive and the false negative errors
#RETURN:
#estimated.probs: estimated marginal, joint and conditional probabilities
"estimate.dag.probs" <-
function( dataset, marginal.probs, joint.probs, parents.pos, error.rates ) {
    #structure where to save the probabilities to be estimated
    estimated.marginal.probs = array(-1, dim=c(nrow(marginal.probs),1));
    estimated.joint.probs = array(-1, dim=c(nrow(marginal.probs),nrow(marginal.probs)));
    estimated.conditional.probs = parents.pos;
    #compute the probability of the AND parents
    parents.probs = array(0,c(nrow(parents.pos),1));
    children.probs = array(0,c(nrow(parents.pos),1));
    last.parent.pos = array(-1,c(nrow(parents.pos),1));
    for (i in 1:nrow(dataset)) {
		for (j in 1:length(parents.pos)) {
			if((length(parents.pos[[j,1]])==1 && parents.pos[[j,1]]==-1)||(sum(dataset[i,parents.pos[[j,1]]])==length(parents.pos[[j,1]]))) {
				parents.probs[j,1] = parents.probs[j,1] + 1;
			}
			if(length(parents.pos[[j,1]])!=1 || parents.pos[[j,1]]!=-1) {
				curr.last.parent = which.min(marginal.probs[parents.pos[[j,1]],1]);
				last.parent.pos[j,1] = curr.last.parent[1];
			}
			if((length(parents.pos[[j,1]])==1 && parents.pos[[j,1]]==-1 && dataset[i,j]==1)||((parents.pos[[j,1]]!=-1) && sum(dataset[i,c(parents.pos[[j,1]],j)])==length(c(parents.pos[[j,1]],j)))) {
				children.probs[j,1] = children.probs[j,1] + 1;
			}
		}
    }
    parents.probs = parents.probs/nrow(dataset);
    children.probs = children.probs/nrow(dataset);
    #estimate the theoretical conditional probabilities given the error rates
    #this estimation is performed by applying the error rates to the marginal and joint probabilities
    theoretical.conditional.probs = array(-1, dim=c(nrow(marginal.probs),1));
    for (i in 1:nrow(theoretical.conditional.probs)) {
        #if the node has a parent, use the error rates to compute the conditional probability
        #if the node has no parent, its conditional probability is not considered
        if(parents.probs[i,1]!=1) {
            #P(i|j)_theoretical = ((P(i,j)_obs-e_p*(P(j)_obs+P(i)_obs)+e_p^2)/(1-e_n-e_p)^2)/((P(j)_obs-e_p)/(1-e_n-e_p))
            theoretical.conditional.probs[i,1] = (children.probs[i,1]-error.rates$error.fp*(parents.probs[i,1]+marginal.probs[i,1])+error.rates$error.fp^2)/((parents.probs[i,1]-error.rates$error.fp)*(1-error.rates$error.fn-error.rates$error.fp));
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
    child.list <- which(parents.probs==1);
    estimated.marginal.probs[child.list,1] = marginal.probs[child.list,1];
    estimated.marginal.probs.with.error = array(-1, dim=c(nrow(marginal.probs),1));
    estimated.marginal.probs.with.error[child.list,1] = estimated.marginal.probs[child.list,1];
    visited = length(child.list);
    #I do not have any constraint for the joint probabilities on the pair of nodes which are the roots of the dag
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
            curr.child <- which(last.parent.pos==node);
            #go through the current children
            for (child in curr.child) {
                #set the marginal probability for this node
                #P(child)_estimate = P(parent)_estimate * P(child|parent)_theoretical
                estimated.marginal.probs[child,1] = estimated.marginal.probs[last.parent.pos[child,1],1]*theoretical.conditional.probs[child,1];
                visited = visited + 1;
                #P(child,parent)_estimare = P(child)_estimate;
                estimated.joint.probs[child,last.parent.pos[child,1]] = estimated.marginal.probs[child,1];
                estimated.joint[child,last.parent.pos[child,1]] = 1;
                estimated.joint.probs[last.parent.pos[child,1],child] = estimated.marginal.probs[child,1];
                estimated.joint[last.parent.pos[child,1],child] = 1;
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
                estimated.joint.probs[i,j] = estimate.dag.joint.probs(i,j,parents.pos,estimated.marginal.probs,theoretical.conditional.probs);
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
    for (i in 1:length(estimated.conditional.probs)) {
		curr.parents.pos = parents.pos[[i,1]];
		for (j in 1:length(parents.pos[[i,1]])) {
			#if the node has no parent, its conditional probability is set to 1
			if(length(curr.parents.pos)==1 && curr.parents.pos==-1) {
				curr.parents.pos = 1;
				break;
			}
			else {
				if(estimated.marginal.probs[curr.parents.pos[j],1]>0) {
					curr.parents.pos[j] = estimated.joint.probs[curr.parents.pos[j],i]/estimated.marginal.probs[curr.parents.pos[j],1];
				}
				else {
					curr.parents.pos[j] = 0;
				}
			}	
        }
        estimated.conditional.probs[[i,1]] = curr.parents.pos;
    }
    #structure to save the results
    estimated.probs = list(marginal.probs=estimated.marginal.probs,joint.probs=estimated.joint.probs,conditional.probs=estimated.conditional.probs);
    return(estimated.probs);
}

#### end of file -- estimate.dag.probs.R

#### estimate.dag.samples.R
####
#### TRONCO: a tool for TRanslational ONCOlogy
####
#### See the files COPYING and LICENSE for copyright and licensing
#### information.


#estimate the probability of observing each sample in the dataset given the reconstructed topology
#INPUT:
#dataset: a valid dataset
#reconstructed.topology: the reconstructed topology
#estimated.marginal.probabilities: estimated marginal probabilities of the events
#estimated.conditional.probabilities: estimated conditional probabilities of the events
#parents.pos: position of the parents of each node
#error.rates: error rates for false positives and false negatives
#RETURN:
#probabilities: probability of each sample
estimate.dag.samples = function(dataset, 
                                reconstructed.topology, 
                                estimated.marginal.probabilities, 
                                estimated.conditional.probabilities, 
                                parents.pos, error.rates) 
{
    #structure where to save the probabilities of the samples
    probabilities = array(-1, c(nrow(dataset), 1))
    
    #compute the position of the latest parent and its conditional probability for each node
    last.parent.pos = array(-1, c(nrow(parents.pos), 1))
    curr.estimated.conditional.probabilities = array(1, c(nrow(estimated.conditional.probabilities), 1))
    
    for (i in 1:length(parents.pos)) {
        if(length(parents.pos[[i, 1]]) != 1 || parents.pos[[i, 1]] != -1) {
            curr.last.parent = which.min(estimated.marginal.probabilities[parents.pos[[i, 1]], 1])
            last.parent.pos[i, 1] = parents.pos[[i, 1]][curr.last.parent[1]]
            curr.estimated.conditional.probabilities[i, 1] = estimated.conditional.probabilities[[i, 1]][curr.last.parent[1]]
        }
    }
    
    #topological properties:
    #1. progression number
    #2. latest parent
    #3. level in the progression
    
    topology.structure = array(0, c(nrow(reconstructed.topology), 3))
    
    #go through the subtrees within the topology
    progression.count = 0
    
    for (i in 1:nrow(reconstructed.topology)) {
        
        #if node i has no parents, it is a root
        if(length(which(reconstructed.topology[, i] == 1)) == 0) {
            progression.count = progression.count + 1
            level = 1
            
            #set the parameters for the root
            topology.structure[i,1] = progression.count
            topology.structure[i,2] = -1
            topology.structure[i,3] = level
            curr.node = i
            
            #go through this progression
            while (length(curr.node) > 0) {
                
                #move to the next level
                level = level + 1
                new.node = vector()
                
                for (j in 1:length(curr.node)) {
                    curr.new.node = which(reconstructed.topology[curr.node[j], ] == 1)
                    
                    if(length(curr.new.node) > 0) {
                        new.node = c(new.node,curr.new.node)
                        
                        for (k in 1:length(curr.new.node)) {
                            
                            #number of the current subprogression
                            topology.structure[curr.new.node[k], 1] = progression.count
                            
                            #parent of the current node
                            if(last.parent.pos[curr.new.node[k], 1] == curr.node[j]) {
                                topology.structure[curr.new.node[k], 2] = curr.node[j]
                            }
                            
                            #level of this node
                            topology.structure[curr.new.node[k], 3] = level
                        }
                    }
                }
                curr.node = new.node
            }
        }
    }

    #go through the dataset and evalutate the probability of each sample
    for (i in 1:nrow(dataset)) {
        sample.probability = 1

        for (j in 1:progression.count) {
            
            #probability of this subprogression (without any knowledge, I set it to 1)
            curr.sample.probability = 1
            
            #entries referring to this subprogression
            curr.entry = which(topology.structure[, 1] == j)
            
            #samples of each element of this subprogression
            curr.sample = dataset[i, curr.entry]
            
            #parents of each element of this subprogression
            curr.parents = topology.structure[curr.entry, 2]
            
            #level of each element of this subprogression
            curr.levels = topology.structure[curr.entry, 3]
            
            #set the probability as the one of the root of this progression
            curr.sample.probability = curr.sample.probability * estimated.marginal.probabilities[curr.entry[which(curr.levels == 1, arr.ind = TRUE)], 1]
            
            #set the maximum level of this subprogression
            max.level = curr.levels[which.max(curr.levels)]

            #if I have at least one event in this sample
            if(length(curr.sample[curr.sample == 1]) > 0) {

                #visit the nodes starting from the lower level
                is.valid = TRUE

                for (k in max.level:1) {
                    curr.level.nodes = which(curr.levels == k, arr.ind=TRUE)

                    #if I'm not on a root
                    if(k > 1) {
                        curr.level.samples = curr.sample[curr.level.nodes]

                        #if I have at least one event at this level
                        if(length(curr.level.samples[curr.level.samples == 1]) > 0) {
                            
                            #I can not have a child without its parent
                            curr.level.parent = curr.parents[curr.level.nodes]
                            
                            for (p in 1:length(curr.level.parent)) {
                                if(dataset[i, curr.level.parent[p]] == 0 && dataset[i, curr.entry[curr.level.nodes[p]]] == 1) {
                                    is.valid = FALSE
                                    break
                                }
                            }
                        }

                        #if the sample is valid
                        if(is.valid == TRUE) {
                            
                            #add the probability of each edge
                            curr.level.parent = curr.parents[curr.level.nodes]
                            
                            for (p in 1:length(curr.level.parent)) {
                                if(dataset[i, curr.level.parent[p]] == 1 && dataset[i,curr.entry[curr.level.nodes[p]]] == 0) {
                                    curr.sample.probability = curr.sample.probability * (1 - curr.estimated.conditional.probabilities[curr.entry[curr.level.nodes[p]], 1])
                                } else if(dataset[i, curr.level.parent[p]]==1 && dataset[i, curr.entry[curr.level.nodes[p]]] == 1) {
                                    curr.sample.probability = curr.sample.probability * curr.estimated.conditional.probabilities[curr.entry[curr.level.nodes[p]], 1]
                                }
                            }
                        }
                    }

                    if(is.valid == FALSE) {
                        curr.sample.probability = 0
                        break
                    }
                }
                if(is.valid == FALSE) {
                    sample.probability = 0
                    break
                }
            } else {
                # ..if this sample has no events for this progression
                curr.sample.probability = 1 - curr.sample.probability
            }

            #update the probability of the topology with the one of this sample
            sample.probability = sample.probability * curr.sample.probability

            if(sample.probability == 0) {
                break
            }
        }
        probabilities[i, 1] = sample.probability
    }

    #correct the estimation by the error rates
    errors.matrix = array(0, c(nrow(probabilities), nrow(dataset)))
    for (i in 1:nrow(probabilities)) {
        for (j in 1:nrow(dataset)) {
            curr.sample.x = as.numeric(dataset[i, ])
            curr.sample.y = as.numeric(dataset[j, ])
            errors.matrix[i, j] = (1 - error.rates$error.fp) ^ ((1 - curr.sample.x) %*% (1 - curr.sample.y)) *
                                error.rates$error.fp ^ ((1 - curr.sample.x) %*% curr.sample.y) *
                                (1 - error.rates$error.fn) ^ (curr.sample.x %*% curr.sample.y) *
                                error.rates$error.fn ^ (curr.sample.x %*% (1 - curr.sample.y))
        }
    }
    probabilities[, 1] = as.numeric(as.vector(probabilities) %*% errors.matrix)
    return(probabilities)
}

#### end of file -- estimate.dag.samples.R
