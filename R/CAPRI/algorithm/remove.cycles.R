#### remove.cycles.R
####
#### TRONCO: a tool for TRanslational ONCOlogy
####
#### See the files COPYING and LICENSE for copyright and licensing
#### information.


# remove any cycle from a given cyclic topology
# INPUT:
# adj.matrix: adjacency matrix of the topology
# weights.temporal.priority: weighted matrix to be used to remove the cycles involving atomic events
# weights.matrix: weighted matrix to be used to remove the cycles involving hypotheses
# not.ordered: list of the nodes to be orderd
# hypotheses: hypotheses to evaluate potential cycles
# RETURN:
# acyclic.topology: structure representing the best acyclic topology
"remove.cycles" <-
function( adj.matrix, weights.temporal.priority, weights.matrix, not.ordered, hypotheses = NA ) {
	
	# evaluate the possible cycles involving atomic events
    if(length(not.ordered)>0) {
    	
    		# consider only the edges that were not ordered by temporal priority
    		curr.edge.pos = 0;
    		for(i in 1:length(not.ordered)) {
    			
    			# consider the events i and j
        		curr.edge = not.ordered[[i]];
        		curr.edge.i = curr.edge[1,1];
        		curr.edge.j = curr.edge[2,1];
        		
        		# check if i and j still create a cycle
        		if(adj.matrix[curr.edge.i,curr.edge.j]==1 && adj.matrix[curr.edge.j,curr.edge.i]==1) {
        			
        			# get the scores of the two edges
        			curr.score.i.j = weights.temporal.priority[curr.edge.i,curr.edge.j];
        			curr.score.j.i = weights.temporal.priority[curr.edge.j,curr.edge.i];
        			
        		# choose an edge based on the score
        		if(curr.score.i.j<curr.score.j.i) {
        			# if i --> j is more confident (lower score) then j --> i
            		cat("Removing edge ",colnames(adj.matrix)[curr.edge.j]," to ",colnames(adj.matrix)[curr.edge.i],"\n");
            		adj.matrix[curr.edge.j,curr.edge.i] = 0;
        		}
        		else {
        			# otherwise
            		cat("Removing edge ",colnames(adj.matrix)[curr.edge.i]," to ",colnames(adj.matrix)[curr.edge.j],"\n");
            		adj.matrix[curr.edge.i,curr.edge.j] = 0;
        		}
        		}
        		
    		}
    		
    	}
    
    # create the structures where to save the weights in increasing order of confidence
    ordered.weights <- vector();
    ordered.edges <- list();
    
    # consider the patterns related the hypotheses
	if(!is.na(hypotheses[1])) {
		
		# evaluate any hypothesis
		res = hypothesis.evaluate.cycles(hypotheses,adj.matrix,as.patterns(hypotheses),weights.matrix);
		
		# matomic provides a map from the atomic events to the hypothesis
		matomic = res$matomic;
		
		# mhypotheses provides a map from the hypotheses to the atomic events
		mhypotheses = res$mhypotheses;
	
		# if I have hypotheses, add the edges to be evaluated during the loop removal
		curr.edge.pos = 0;
		for(i in 1:nrow(adj.matrix)) {
			for(j in 1:nrow(adj.matrix)) {
				if(adj.matrix[i,j]==1) {
					ordered.weights = rbind(ordered.weights,weights.matrix[i,j]);
					curr.edge.pos = curr.edge.pos + 1;
            		new.edge <- array(0, c(2,1));
            		new.edge[1,1] = i;
            		new.edge[2,1] = j;
            		ordered.edges[curr.edge.pos] = list(new.edge);
				}
			}
		}
	
		# sort the edges in increasing order of confidence (i.e. the edges with lower pvalue are the most confident)
    	ordered.edges = ordered.edges[sort(unlist(ordered.weights),decreasing=TRUE,index.return=TRUE)$ix];
	
	}
    
    # visit the ordered edges and remove the ones that are causing any cycle
    if(length(ordered.edges)>0) {
    	
    		total.edges = length(which(adj.matrix == 1))
    		removed = 0
        for(i in 1:length(ordered.edges)) {
        	
            # consider the edge i-->j
            curr.edge = ordered.edges[[i]];
            curr.edge.i = curr.edge[1,1];
            curr.edge.j = curr.edge[2,1];
            
            # search for loops between curr.edge.i and curr.edge.j
            is.path = find.path(adj.matrix,curr.edge.j,curr.edge.i,matomic,mhypotheses,vector());
            
            # if there is a path between the two nodes, remove edge i --> j
            if(is.path==1) {
				removed = removed + 1
            		# cat("Removing edge ",colnames(adj.matrix)[curr.edge.i]," to ",colnames(adj.matrix)[curr.edge.j],"\n");
            		adj.matrix[curr.edge.i,curr.edge.j] = 0;
            }
            
        }
        
        cat(paste0('\tRemoved ', removed, ' edges out of ', total.edges ,' (', round(100 * removed/total.edges, 0),'%)\n'))
    }
    
    # save the results and return them
    acyclic.topology = list(adj.matrix=adj.matrix);
    return(acyclic.topology);
    
}

#### end of file -- remove.cycles.R
