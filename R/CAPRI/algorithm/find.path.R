#### find.path.R
####
#### TRONCO: a tool for TRanslational ONCOlogy
####
#### See the files COPYING and LICENSE for copyright and licensing
#### information.


# search for a path between the two given nodes
# INPUT:
# adj.matrix: adjacency matrix of the topology
# first.node: first node of the path
# last.node: last node of the path
# hypotheses: all the hypotheses
# hatomic: map of the hypotheses and their atomic events
# visited.nodes: list of all the visited nodes
# RETURN:
# is.path: 0 if there is not a path from first.node to last.node, 1 if there is
"find.path" <-
function(adj.matrix, first.node, last.node, hypotheses = NA, hatomic = NA, visited.nodes) {
	
	# suppose that there is no path from first.node to last.node
    is.path = 0;
    
    # if hypotheses, expand first.node and last.node to their atomic events
    curr.first.node = first.node;
    if(exists(toString(first.node), envir = hatomic)) {
    	curr.first.node = hatomic[[toString(first.node)]];
    }
    curr.last.node = last.node;
    if(exists(toString(last.node), envir = hatomic)) {
    	curr.last.node = hatomic[[toString(last.node)]];
    }
    
    # recursion base case: first.node==last.node
    if(any(curr.first.node%in%curr.last.node)) {
        is.path = 1;
    }
    
    # start the recursive search
    else {
    	
        # find the possible paths starting from first.node
        # [first.node,] are the nodes directly caused by first.node
        possible.paths = which(adj.matrix[first.node,]==1);
        
        # expand the possible paths considering the atomic events of each hypothesis
        if(length(possible.paths)>0 && length(hatomic)>0) {
        		curr.atomic.pool = vector();
        		for(i in 1:length(possible.paths)) {
        			# get any atomic element of the current hypothesis
        			if(exists(toString(possible.paths[i]), envir = hatomic)) {
        				curr.atomic.pool = unique(c(curr.atomic.pool,hatomic[[toString(possible.paths[i])]]));
        			}
        		}
        		possible.paths = unique(c(possible.paths,curr.atomic.pool));
        }
        
        # expand the possible paths considering the hypotheses of each atomic event
        if(length(possible.paths)>0 && !is.na(hypotheses[1])) {
        		curr.hypotheses.pool = vector();
        		for(i in 1:length(possible.paths)) {
        			# get any hypothesis of the current atomic element
        			new.hypotheses.pool = events.pattern(hypotheses,colnames(adj.matrix)[possible.paths[i]]);
        			if(!is.na(new.hypotheses.pool) && length(new.hypotheses.pool)>0) {
        				for(j in 1:length(new.hypotheses.pool)) {
        					# add the hypothesis
        					curr.hypotheses.pool = unique(c(curr.hypotheses.pool,which(colnames(adj.matrix)==new.hypotheses.pool[j])));
        				# add the atoms in the hypothesis
        				if(exists(toString(which(colnames(adj.matrix)==new.hypotheses.pool[j])), envir = hatomic)) {
        					curr.hypotheses.pool = unique(c(curr.hypotheses.pool,hatomic[[toString(which(colnames(adj.matrix)==new.hypotheses.pool[j]))]]));
        				}
        				}
        			}
        		}
        		possible.paths = unique(c(possible.paths,curr.hypotheses.pool));
        }
        
        # recursive case: there are paths starting from first.node
        if(length(possible.paths)>0) {
            curr.path = 0;
            while(curr.path<length(possible.paths) && is.path==0) {
                curr.path = curr.path + 1;
                # check if I'm in a cycle
                if(length(visited.nodes)>0 && length(which(visited.nodes==possible.paths[curr.path]))>0) {
                    break;
                }
                # if I'm not, go on
                else {
                		visited.nodes = rbind(visited.nodes,possible.paths[curr.path]);
                    if(possible.paths[curr.path]==last.node) {
                        is.path = 1;
                    }
                    else {
                        find.path(adj.matrix,possible.paths[curr.path],last.node,hypotheses,hatomic,visited.nodes);
                    }
                }
            }
        }
        
    }
    
    # return the results
    return(is.path);
    
}

#### end of file -- find.path.R
