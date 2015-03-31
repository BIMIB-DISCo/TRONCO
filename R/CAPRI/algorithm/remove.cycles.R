#### remove.cycles.R
####
#### TRONCO: a tool for TRanslational ONCOlogy
####
#### See the files COPYING and LICENSE for copyright and licensing
#### information.


# remove any cycle from a given cyclic topology
# INPUT:
# adj.matrix: adjacency matrix of the topology
# weights.matrix: weighted matrix to be used to remove the cycles
# not.ordered: list of the nodes to be orderd
# hypotheses: hypotheses to evaluate potential cycles
# RETURN:
# acyclic.topology: structure representing the best acyclic topology
"remove.cycles" <-
function( adj.matrix, weights.matrix, not.ordered, hypotheses = NA ) {
	
	# create the structures where to save the weights in increasing order of confidence
    ordered.weights <- vector();
    ordered.edges <- list();
    
    # evaluate the possible cycles involving atomic events
    if(length(not.ordered)>0) {
    	
    		# consider only the edges that were not ordered by temporal priority
    		curr.edge.pos = 0;
    		for(i in 1:length(not.ordered)) {
    			
    			# consider the events i and j
        		curr.edge = not.ordered[[i]];
        		curr.edge.i = curr.edge[1,1];
        		curr.edge.j = curr.edge[2,1];
        		
        		# check if i --> j is a valid edge
        		if(adj.matrix[curr.edge.i,curr.edge.j]==1) {
            		ordered.weights = rbind(ordered.weights,weights.matrix[curr.edge.i,curr.edge.j]);
            		curr.edge.pos = curr.edge.pos + 1;
            		new.edge <- array(0, c(3,1));
            		new.edge[1,1] = curr.edge.i;
            		new.edge[2,1] = curr.edge.j;
            		ordered.edges[curr.edge.pos] = list(new.edge);
        		}
        		
        		# check if j --> i is a valid edge
        		if(adj.matrix[curr.edge.j,curr.edge.i]==1) {
            		ordered.weights = rbind(ordered.weights,weights.matrix[curr.edge.j,curr.edge.i]);
            		curr.edge.pos = curr.edge.pos + 1;
            		new.edge <- array(0, c(3,1));
            		new.edge[1,1] = curr.edge.j;
            		new.edge[2,1] = curr.edge.i;
            		ordered.edges[curr.edge.pos] = list(new.edge);
        		}
        		
    		}
    		
    	}
    
    # consider the patterns related the hypotheses
	if(!is.na(hypotheses[1])) {
		
		# evaluate any hypothesis
		res = hypothesis.evaluate.cycles(hypotheses,adj.matrix,as.patterns(hypotheses),weights.matrix);
		
		# of each pattern, compute the involving edges and weights
		ordered.edges = c(ordered.edges,res$ordered.edges);
		ordered.weights = c(ordered.weights,res$ordered.weights);
		
		# hatomic provides a map of the atomic events of any hypothesis
		hatomic = res$hatomic;
	
	}
    
    # sort the edges in increasing order of confidence (i.e. the edges with lower pvalue are the most confident)
    ordered.edges = ordered.edges[sort(unlist(ordered.weights),decreasing=TRUE,index.return=TRUE)$ix];
    
    # visit the ordered edges and remove the ones that are causing cycles, if any
    if(length(ordered.edges)>0) {
        for(i in 1:length(ordered.edges)) {
        	
            # consider the edge i-->j
            curr.edge = ordered.edges[[i]];
            curr.edge.i = curr.edge[1,1];
            curr.edge.j = curr.edge[2,1];
            
            # if this edge does not involve any hypothesis
            if(curr.edge[3,1]==0) {
            		is.path = find.path(adj.matrix,curr.edge.j,curr.edge.i,hypotheses,hatomic,vector());
            }
            
            # otherwise, if this edge does involve an hypothesis
            else {
            		
            		# get the atomic elements of the current hypothesis
            		atomic.pool = hatomic[[toString(curr.edge[curr.edge[3,1],1])]];
            		
            		# navigate the atomic elements searching for loops 
            		for(j in 1:length(atomic.pool)) {
            			
            			# search for any loop on the incoming at outgoing edges involving any hypothesis
            			is.path = 0;
            			if(curr.edge[3,1]==1) {
            				is.path = find.path(adj.matrix,curr.edge.j,atomic.pool[j],hypotheses,hatomic,vector());
            			}
            			else if(curr.edge[3,1]==2) {
            				is.path = find.path(adj.matrix,atomic.pool[j],curr.edge.i,hypotheses,hatomic,vector());
            			}
            			
            			# if I find the first loop, stop searching
            			if(is.path==1) {
            				break;
            			}
            			
            		}
            		
            }
            
            # if there is a path between the two nodes, remove edge i --> j
            if(is.path==1) {
                adj.matrix[curr.edge.i,curr.edge.j] = 0;
            }
            
        }
    }
    
    # save the results and return them
    acyclic.topology = list(adj.matrix=adj.matrix);
    return(acyclic.topology);
    
}

#### end of file -- remove.cycles.R
