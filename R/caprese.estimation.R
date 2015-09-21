##################################################################################
#                                                                                #
# TRONCO: a tool for TRanslational ONCOlogy                                      #
#                                                                                #
##################################################################################
# Copyright (c) 2015, Marco Antoniotti, Giulio Caravagna, Luca De Sano,          #
# Alex Graudenzi, Ilya Korsunsky, Mattia Longoni, Loes Olde Loohuis,             #
# Giancarlo Mauri, Bud Mishra and Daniele Ramazzotti.                            #
#                                                                                #
# All rights reserved. This program and the accompanying materials               #
# are made available under the terms of the GNU GPL v3.0                         #
# which accompanies this distribution                                            #
#                                                                                #
##################################################################################


# estimate the error rates by "L-BFGS-B" optimization in terms of L2-error
# @title estimate.tree.error.rates
# @param marginal.probs marginal probabilities
# @param joint.probs joint probabilities
# @param parents.pos which event is the parent? 0 if none, a number otherwise
# @return estimated.error.rates: estimated probabilities, false positive and false negative error rates
estimate.tree.error.rates = function(marginal.probs,
                                     joint.probs,
                                     parents.pos) 
{
    
    # function to be optimized by "L-BFGS-B" optimization in terms of L2-error
    f.estimation <- function(errors) {
        # set the current error rates with the starting point of the optimization being (e_pos,e_neg) = 0
        error.rates = list(error.fp=errors[1], error.fn=errors[2])
        # estimate the observed probabilities given the error rates
        estimated.probs = estimate.tree.probs(marginal.probs, joint.probs, parents.pos,error.rates)
        # evaluate the goodness of the estimatione by L2-error on the estimated marginal and joint probabilities
        error.estimation = sum((marginal.probs - estimated.probs$marginal.probs)^2) + sum((joint.probs-estimated.probs$joint.probs)^2)
        return(error.estimation);
    }
    
    # the estimation is performed as in Byrd et al (1995)
    # this method allows for box constraints, i.e., each variable can be given a lower and/or upper bound
    estimated.error.rates = optim(c(0.00, 0.00),
                                  f.estimation,
                                  method="L-BFGS-B",
                                  lower=c(0.00, 0.00),
                                  upper=c(0.49, 0.49))$par
    
    # structure to save the results
    estimated.error.rates = list(error.fp=estimated.error.rates[1], error.fn=estimated.error.rates[2])
    return(estimated.error.rates)
    
}

# estimate the theoretical joint probability of two given nodes given the reconstructed topology
# @title estimate.tree.joint.probs
# @param first.node first node
# @param second.node second node
# @param parents.pos which event is the parent? -1 if none, a number otherwise
# @param marginal.probs marginal probabilities
# @param conditional.probs conditional probabilities
# @return estimated.tree.joint.probs: estimated theoretical joint probability
estimate.tree.joint.probs = function(first.node,
                                     second.node,
                                     parents.pos,
                                     marginal.probs,
                                     conditional.probs) 
{
    
    # if the two nodes are roots
    if(parents.pos[first.node] == -1 && parents.pos[second.node] == -1) {
        estimated.tree.joint.probs = marginal.probs[first.node, 1] * marginal.probs[second.node, 1]
    } else if(first.node == second.node) {
        # if the two nodes are the same node
        estimated.tree.joint.probs = marginal.probs[first.node, 1]
    } else {
        # otherwise
        # go through the parents starting from the two nodes to find if they are directly connected
        # are the two nodes in the same path?
        is.path = 0;
        # check if first.node is an ancestor of second.node
        curr.first = first.node
        curr.second = second.node
        while(parents.pos[curr.second] != -1) {
            if(curr.first == curr.second) {
                is.path = 1
                is.child = second.node
                break
            }
            curr.second = parents.pos[curr.second]
        }
        if(is.path == 0) {
            # check if second.node is an ancestor of first.node
            curr.first = first.node;
            curr.second = second.node;
            while(parents.pos[curr.first] != -1) {
                if(curr.first == curr.second) {
                    is.path = 1
                    is.child = first.node
                    break
                }
                curr.first = parents.pos[curr.first]
            }
        }
        # check if the two nodes are at least connected
        # are the two nodes connected?
        is.connected = 0
        if(is.path == 0) {
            curr.first = first.node
            curr.second = second.node
            while(parents.pos[curr.first] != -1) {
                while(parents.pos[curr.second] != -1) {
                    if(curr.first == curr.second) {
                        is.connected = 1
                        is.ancestor = curr.first
                        break
                    } else {
                        curr.second = parents.pos[curr.second]
                    }
                }
                if(is.connected == 0) {
                    curr.first = parents.pos[curr.first]
                    curr.second = second.node
                } else {
                    break;
                }
            }
        }
        # now I can set the joint probabilities
        # in this case the two nodes are directly connected
        # P(child,parent)_estimate = P(child);
        if(is.path == 1) {
            estimated.tree.joint.probs = marginal.probs[is.child, 1]
        } else if(is.connected == 1) {
            # in this case the two nodes are indirectly connected
            # P(i,j)_estimate = P(ancestor)_estimate * P_PATH(ancestor->first.node)_estimate * P_PATH(ancestor->second.node)_estimate
            # P(ancestor)_estimate
            estimated.tree.joint.probs = marginal.probs[is.ancestor, 1]
            # P_PATH(ancestor->first.node)_estimate
            first.path = 1
            curr.first = first.node
            while(parents.pos[curr.first] != is.ancestor) {
                first.path = first.path * conditional.probs[curr.first, 1]
                curr.first = parents.pos[curr.first]
            }
            second.path = 1
            curr.second = second.node
            while(parents.pos[curr.second] != is.ancestor) {
                second.path = second.path * conditional.probs[curr.second,1]
                curr.second = parents.pos[curr.second]
            }
            estimated.tree.joint.probs = estimated.tree.joint.probs * first.path * second.path
        } else {
            # in this case the two nodes are not connected
            # P(i,j)_estimate = P(i)_estimate * P(j)_estimate
            estimated.tree.joint.probs = marginal.probs[first.node, 1] * marginal.probs[second.node, 1]
        }
    }
    return(estimated.tree.joint.probs)
    
}

# estimate the marginal, joint and conditional probabilities given the reconstructed topology and the error rates
# @title estimate.tree.probs
# @param marginal.probs observed marginal probabilities
# @param joint.probs observed joint probabilities
# @param parents.pos position of the parents in the list of nodes
# @param error.rates rates for the false positive and the false negative errors
# @return estimated.probs estimated marginal, joint and conditional probabilities
estimate.tree.probs = function(marginal.probs,
                               joint.probs,
                               parents.pos,
                               error.rates) {
    
    # structure where to save the probabilities to be estimated
    estimated.marginal.probs = array(-1, dim=c(nrow(marginal.probs), 1))
    estimated.joint.probs = array(-1, dim=c(nrow(marginal.probs), nrow(marginal.probs)))
    estimated.conditional.probs = array(-1, dim=c(nrow(marginal.probs), 1))
    
    # estimate the theoretical conditional probabilities given the error rates
    # this estimation is performed by applying the error rates to the marginal and joint probabilities
    theoretical.conditional.probs = array(-1, dim=c(nrow(marginal.probs),1))
    for (i in 1:nrow(theoretical.conditional.probs)) {
        # if the node has a parent, use the error rates to compute the conditional probability
        # if the node has no parent, its conditional probability is not considered
        if(parents.pos[i,1] != -1) {
            # P(i|j)_theoretical = ((P(i,j)_obs-e_p*(P(j)_obs+P(i)_obs)+e_p^2)/(1-e_n-e_p)^2)/((P(j)_obs-e_p)/(1-e_n-e_p))
            theoretical.conditional.probs[i,1] = (joint.probs[i, parents.pos[i,1]] - error.rates$error.fp * (marginal.probs[parents.pos[i,1],1]+marginal.probs[i,1])+error.rates$error.fp^2)/((marginal.probs[parents.pos[i,1],1]-error.rates$error.fp)*(1-error.rates$error.fn-error.rates$error.fp));
            if(theoretical.conditional.probs[i,1] < 0 || theoretical.conditional.probs[i,1]>1) {
                # invalid theoretical conditional probability
                if(theoretical.conditional.probs[i,1] < 0) {
                    theoretical.conditional.probs[i,1] = 0
                } else {
                    theoretical.conditional.probs[i,1] = 1
                }
            }
        }
    }
    
    # estimate the marginal observed probabilities
    # this estimation is performed by applying the topological constraints on the probabilities and then the error rates
    # I do not have any constraint on the nodes without a parent
    child.list = which(parents.pos == -1)
    estimated.marginal.probs[child.list, 1] = marginal.probs[child.list, 1]
    estimated.marginal.probs.with.error = array(-1, dim=c(nrow(marginal.probs), 1))
    estimated.marginal.probs.with.error[child.list, 1] = estimated.marginal.probs[child.list, 1]
    visited = length(child.list)
    # I do not have any constraint for the joint probabilities on the pair of nodes which are the roots of the tree/forest
    estimated.joint = array(0, dim=c(nrow(marginal.probs), nrow(marginal.probs)))
    for (i in child.list) {
        for (j in child.list) {
            if(i!=j) {
                estimated.joint.probs[i, j] = joint.probs[i, j]
                estimated.joint[i,j] = -1
            }
        }
    }
    # visit the nodes with a parent in topological order
    while (visited < nrow(estimated.marginal.probs)) {
        # set the new child list
        new.child = vector()
        # go through the current parents
        for (node in child.list) {
            # set the new children
            curr.child <- which(parents.pos == node)
            # go through the current children
            for (child in curr.child) {
                # set the marginal probability for this node
                # P(child)_estimate = P(parent)_estimate * P(child|parent)_theoretical
                estimated.marginal.probs[child,1] = estimated.marginal.probs[parents.pos[child, 1], 1] * theoretical.conditional.probs[child,1];
                visited = visited + 1
                # P(child,parent)_estimare = P(child)_estimate;
                estimated.joint.probs[child,parents.pos[child, 1]] = estimated.marginal.probs[child, 1]
                estimated.joint[child,parents.pos[child, 1]] = 1
                estimated.joint.probs[parents.pos[child, 1], child] = estimated.marginal.probs[child, 1]
                estimated.joint[parents.pos[child, 1], child] = 1
                # apply the error rates to the marginal probabilities
                # P(i)_obs_estimate = P(i)_estimate*(1-e_n) + P(not i)_estimate*e_p
                estimated.marginal.probs.with.error[child, 1] = error.rates$error.fp + (1 - error.rates$error.fn-error.rates$error.fp) * estimated.marginal.probs[child, 1]
                if(estimated.marginal.probs.with.error[child, 1] < 0 || estimated.marginal.probs.with.error[child, 1] > 1) {
                    # invalid estimated observed probability
                    if(estimated.marginal.probs.with.error[child, 1] < 0) {
                        estimated.marginal.probs.with.error[child, 1] = 0
                    } else {
                        estimated.marginal.probs.with.error[child, 1] = 1
                    }
                }
            }
            new.child = c(new.child,curr.child)
        }
        # set the next child list
        child.list = new.child
    }
    diag(estimated.joint.probs) = estimated.marginal.probs
    diag(estimated.joint) = -1
    
    # given the estimated observed probabilities, I can now also estimate the joint probabilities by applying the topological constraints and then the error rates
    for (i in 1:nrow(estimated.joint.probs)) {
        for (j in i:nrow(estimated.joint.probs)) {
            # if I still need to estimate this joint probability
            if(estimated.joint[i,j] == 0) {
                estimated.joint.probs[i, j] = estimate.tree.joint.probs(i, 
                                                                        j, 
                                                                        parents.pos, 
                                                                        estimated.marginal.probs,
                                                                        theoretical.conditional.probs)
                estimated.joint[i, j] = 1
            }
            # now I can apply the error rates to estimate the observed joint probabilities
            if(estimated.joint[i,j] == 1) {
                # P(i,j)_obs_estimate = P(i,j)_estimate*(1-e_n)^2+P(not i,j)_estimate*e_p*(1-e_n)+P(i,not j)_estimate*(1-e_n)*e_p+P(not i,not j)_estimate*e_p^2;
                estimated.joint.probs[i,j] = estimated.joint.probs[i, j] * ((1-error.rates$error.fn - error.rates$error.fp)^2) + error.rates$error.fp * (estimated.marginal.probs[i, 1] + estimated.marginal.probs[j, 1]) - error.rates$error.fp^2
                # invalid estimated joint probability
                if(estimated.joint.probs[i,j] < 0 || estimated.joint.probs[i, j] > min(estimated.marginal.probs.with.error[i, 1], estimated.marginal.probs.with.error[j, 1])) {
                    if(estimated.joint.probs[i,j] < 0) {
                        estimated.joint.probs[i,j] = 0
                    } else {
                        estimated.joint.probs[i, j] = min(estimated.marginal.probs.with.error[i, 1], estimated.marginal.probs.with.error[j, 1])
                    }
                }
                estimated.joint.probs[j, i] = estimated.joint.probs[i,j]
            }
        }
    }
    
    # save the estimated probabilities
    estimated.marginal.probs = estimated.marginal.probs.with.error
    # given the estimated observed and joint probabilities, I can finally compute the conditional probabilities
    # P(child|parent)_obs_estimate = P(parent,child)_obs_estimate/P(parent)_obs_estimate
    for (i in 1:nrow(estimated.conditional.probs)) {
        if(parents.pos[i, 1] != -1) {
            if(estimated.marginal.probs[parents.pos[i, 1], 1] > 0) {
                estimated.conditional.probs[i, 1] = estimated.joint.probs[parents.pos[i, 1], i] / estimated.marginal.probs[parents.pos[i, 1], 1]
            } else {
                estimated.conditional.probs[i, 1] = 0
            }
        } else {
            # if the node has no parent, its conditional probability is set to 1
            estimated.conditional.probs[i, 1] = 1
        }
    }
    
    # structure to save the results
    estimated.probs = list(marginal.probs=estimated.marginal.probs,
                           joint.probs=estimated.joint.probs,
                           conditional.probs=estimated.conditional.probs)
    return(estimated.probs);
    
}

# estimate the probability of observing each sample in the dataset given the reconstructed topology
# @title estimate.tree.samples
# @param dataset a valid dataset
# @param reconstructed.topology the reconstructed topology
# @param estimated.marginal.probabilities estimated marginal probabilities of the events
# @param estimated.conditional.probabilities estimated conditional probabilities of the events
# @param error.rates error rates for false positives and false negatives
# @return probabilities: probability of each sample
estimate.tree.samples = function(dataset,
                                 reconstructed.topology, 
                                 estimated.marginal.probabilities, 
                                 estimated.conditional.probabilities, 
                                 error.rates) {
    
    # structure where to save the probabilities of the samples
    probabilities = array(-1, c(nrow(dataset), 1))
    # topological properties:
    # 1. tree number
    # 2. parent
    # 3. level in the tree
    topology.structure = array(0, c(nrow(reconstructed.topology), 3))
    
    # go through the subtrees within the topology of four 
    tree.count = 0
    for (i in 1:nrow(reconstructed.topology)) {
        # if node i has no parents, it is a root
        if(length(which(reconstructed.topology[,i] == 1)) == 0) {
            tree.count = tree.count + 1
            level = 1
            # set the parameters for the root
            topology.structure[i,1] = tree.count
            topology.structure[i,2] = -1
            topology.structure[i,3] = level
            curr.node = i
            # go through this tree
            while (length(curr.node) > 0) {
                # move to the next level
                level = level + 1
                new.node = vector()
                for (j in 1:length(curr.node)) {
                    curr.new.node = which(reconstructed.topology[curr.node[j], ] == 1)
                    if(length(curr.new.node) > 0) {
                        new.node = c(new.node,curr.new.node)
                        for (k in 1:length(curr.new.node)) {
                            # number of the current subtree
                            topology.structure[curr.new.node[k], 1] = tree.count
                            # parent of the current node
                            topology.structure[curr.new.node[k], 2] = curr.node[j]
                            # level of this node
                            topology.structure[curr.new.node[k], 3] = level
                        }
                    }
                }
                curr.node = new.node
            }
        }
    }
    
    # go through the dataset and evalutate the probability of each sample
    for (i in 1:nrow(dataset)) {
        sample.probability = 1
        for (j in 1:tree.count) {
            # probability of this subtree (without any knowledge, I set it to 1)
            curr.sample.probability = 1
            # entries referring to this subtree
            curr.entry = which(topology.structure[, 1] == j)
            # samples of each element of this subtree
            curr.sample = dataset[i,curr.entry]
            # parents of each element of this subtree
            curr.parents = topology.structure[curr.entry, 2]
            # level of each element of this subtree
            curr.levels = topology.structure[curr.entry, 3]
            # set the probability as the one of the root of this tree
            curr.sample.probability = curr.sample.probability * estimated.marginal.probabilities[curr.entry[which(curr.levels == 1, arr.ind=TRUE)], 1]
            # set the maximum level of this subtree
            max.level = curr.levels[which.max(curr.levels)]
            # if I have at least one event in this sample
            if(length(curr.sample[curr.sample == 1]) > 0) {
                # visit the nodes starting from the lower level
                is.valid = TRUE
                for (k in max.level:1) {
                    curr.level.nodes = which(curr.levels == k, arr.ind=TRUE)
                    # if I'm not on a root
                    if(k > 1) {
                        curr.level.samples = curr.sample[curr.level.nodes]
                        # if I have at least one event at this level
                        if(length(curr.level.samples[curr.level.samples == 1]) > 0) {
                            # I can not have a child without its parent
                            curr.level.parent = curr.parents[curr.level.nodes]
                            for (p in 1:length(curr.level.parent)) {
                                if(dataset[i,curr.level.parent[p]] == 0 && dataset[i, curr.entry[curr.level.nodes[p]]] == 1) {
                                    is.valid = FALSE
                                    break
                                }
                            }
                        }
                        # if the sample is valid
                        if(is.valid == TRUE) {
                            # add the probability of each edge
                            curr.level.parent = curr.parents[curr.level.nodes]
                            for (p in 1:length(curr.level.parent)) {
                                if(dataset[i, curr.level.parent[p]] == 1 && dataset[i,curr.entry[curr.level.nodes[p]]] == 0) {
                                    curr.sample.probability = curr.sample.probability * (1 - estimated.conditional.probabilities[curr.entry[curr.level.nodes[p]], 1])
                                }
                                else if(dataset[i,curr.level.parent[p]]==1 && dataset[i, curr.entry[curr.level.nodes[p]]] == 1) {
                                    curr.sample.probability = curr.sample.probability * estimated.conditional.probabilities[curr.entry[curr.level.nodes[p]], 1]
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
            # if this sample has no events for this tree
                curr.sample.probability = 1 - curr.sample.probability
            }
            # update the probability of the topology with the one of this sample
            sample.probability = sample.probability * curr.sample.probability
            if(sample.probability == 0) {
                break
            }
        }
        probabilities[i, 1] = sample.probability;
    }
    
    # correct the estimation by the error rates
    errors.matrix <- array(0, c(nrow(probabilities), nrow(dataset)))
    for (i in 1:nrow(probabilities)) {
        for (j in 1:nrow(dataset)) {
            curr.sample.x = as.numeric(dataset[i, ])
            curr.sample.y = as.numeric(dataset[j, ])
            errors.matrix[i, j] = (1 - error.rates$error.fp)^((1 - curr.sample.x) %*% (1 - curr.sample.y)) * error.rates$error.fp ^ ((1 - curr.sample.x) %*% curr.sample.y) * (1 - error.rates$error.fn)^(curr.sample.x %*% curr.sample.y) * error.rates$error.fn^(curr.sample.x %*% (1 - curr.sample.y))
        }
    }
    
    probabilities[, 1] = as.numeric(as.vector(probabilities) %*% errors.matrix)
    return(probabilities)
    
}
