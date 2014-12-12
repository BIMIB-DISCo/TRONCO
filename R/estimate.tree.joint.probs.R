##################################################################################
#                                                                                #
# TRONCO: a tool for TRanslational ONCOlogy                                      #
#                                                                                #
##################################################################################
# Copyright (c) 2014, Marco Antoniotti, Giulio Caravagna, Alex Graudenzi,        #
# Ilya Korsunsky, Mattia Longoni, Loes Olde Loohuis, Giancarlo Mauri, Bud Mishra #
# and Daniele Ramazzotti.                                                        #
#                                                                                #
# All rights reserved. This program and the accompanying materials               #
# are made available under the terms of the Eclipse Public License v1.0          #
# which accompanies this distribution, and is available at                       #
# http://www.eclipse.org/legal/epl-v10.html and in the include COPYING file      #
#                                                                                #
# Initial contributors:                                                          #
# Giulio Caravagna, Alex Graudenzi, Mattia Longoni and Daniele Ramazzotti.       #
##################################################################################

#estimate the theoretical joint probability of two given nodes given the reconstructed topology
#INPUT:
#first.node: first node
#second.node: second node
#parents.pos: which event is the parent? 0 if none, a number otherwise
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
                is.child = curr.second;
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
                    is.child = curr.first;
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
            #P(ancestor)_estimat
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
