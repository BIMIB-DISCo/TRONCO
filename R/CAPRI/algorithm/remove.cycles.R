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
            		new.edge <- array(0, c(2,1));
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
		
		# matomic provides a map from the atomic events to the hypothesis
		matomic = res$matomic;
		
		# mhypotheses provides a map from the hypotheses to the atomic events
		mhypotheses = res$mhypotheses;
	
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
            
            # search for loops between curr.edge.i and curr.edge.j
            is.path = find.path(adj.matrix,curr.edge.j,curr.edge.i,matomic,mhypotheses,vector());
            
            # if there is a path between the two nodes, remove edge i --> j
            if(is.path==1) {
            		cat("Removing edge ",colnames(adj.matrix)[curr.edge.i]," to ",colnames(adj.matrix)[curr.edge.j],"\n");
            		adj.matrix[curr.edge.i,curr.edge.j] = 0;
            }
            
        }
        
    }
    
    # save the results and return them
    acyclic.topology = list(adj.matrix=adj.matrix);
    return(acyclic.topology);
    
}

#### end of file -- remove.cycles.R
