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
# visited.nodes: list of all the visited nodes
# RETURN:
# is.path: 0 if there is not a path from first.node to last.node, 1 if there is
"find.path" <-
function(adj.matrix, first.node, last.node, visited.nodes) {
	
    # suppose that there is no path from first.node to last.node
    is.path = 0;
    
    # recursion base case: first.node==last.node
    if(first.node==last.node) {
        is.path = 1;
    }
    # start the recursive search
    else {
        # find the possible paths starting from first.node
        # [first.node,] are the nodes directly caused by first.node
        possible.paths = which(adj.matrix[first.node,]==1);
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
                        find.path(adj.matrix,possible.paths[curr.path],last.node,visited.nodes);
                    }
                }
            }
        }
    }
    
    # return the results
    return(is.path);
    
}

#### end of file -- find.path.R
