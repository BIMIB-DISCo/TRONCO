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
# matomic: map of the hypotheses to their atomic events
# mhypotheses: map of the atomic events to their hypotheses
# visited.nodes: list of all the visited nodes
# RETURN:
# is.path: 0 if there is not a path from first.node to last.node, 1 if there is
"find.path" <-
function( adj.matrix, first.node, last.node, matomic = NA, mhypotheses = NA, visited.nodes ) {
	
	# suppose that there is no path from first.node to last.node
    is.path = 0;
    
    # if first.node is an hypothesis, expand it to its atomic events
    curr.first.node = vector();
    if(exists(toString(first.node), envir = matomic)) {
    		curr.first.node = matomic[[toString(first.node)]];
    }
    
    # otherwise, if first.node is an atomic event, expand it to its hypotheses
    else if(exists(toString(first.node), envir = mhypotheses)) {
		curr.first.node = mhypotheses[[toString(first.node)]];
    }
    
    # append to first.node any alias
    first.node = c(first.node,curr.first.node);
    
    # if last.node is an hypothesis, expand it to its atomic events
    curr.last.node = vector();
    if(exists(toString(last.node), envir = matomic)) {
    		curr.last.node = matomic[[toString(last.node)]];
    }
    
    # otherwise, if last.node is an atomic event, expand it to its hypotheses
    else if(exists(toString(last.node), envir = mhypotheses)) {
		curr.last.node = mhypotheses[[toString(last.node)]];
    }
    
    # append to curr.last.node any alias of last.node
    curr.last.node = c(last.node,curr.last.node);
    
    # recursion base case: first.node==last.node
    if(any(first.node%in%curr.last.node)) {
        is.path = 1;
    }
    
    # start the recursive search
    else {
    		
		# find the possible paths starting from first.node
		# [first.node,] are the nodes directly caused by first.node
    		possible.paths =  vector();
    		if(length(first.node)>0) {
    			for (i in 1:length(first.node)) {
    				possible.paths = unique(c(possible.paths,which(adj.matrix[first.node[i],]==1)));
    			}
    		}
	
		# recursive case: there are paths starting from first.node
		if(length(possible.paths)>0) {
		
			curr.path = 0;
			while(curr.path<length(possible.paths) && is.path==0) {
				curr.path = curr.path + 1;
				# check if I'm in a cycle
				if(length(visited.nodes)>0 && any(which(visited.nodes==possible.paths[curr.path]))) {
					break;
				}
				# if I'm not, go on
				else {
					visited.nodes = rbind(visited.nodes,possible.paths[curr.path]);
					
					if(any(possible.paths[curr.path]==curr.last.node)) {
						is.path = 1;
					}
					else {
						is.path = find.path(adj.matrix,possible.paths[curr.path],last.node,matomic,mhypotheses,visited.nodes);
					}
				}
			}
            
		}
        
    }
    
    # return the results
    return(is.path);
    
}

#### end of file -- find.path.R
