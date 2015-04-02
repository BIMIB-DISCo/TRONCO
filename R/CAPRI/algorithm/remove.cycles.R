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
	
	total.edges = length(which(adj.matrix == 1))
    removed = 0
	
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
	        			removed = removed + 1
	            		# cat("Removing edge ",colnames(adj.matrix)[curr.edge.j]," to ",colnames(adj.matrix)[curr.edge.i],"\n");
	            		adj.matrix[curr.edge.j,curr.edge.i] = 0;
	        		}
	        		else {
	        			# otherwise
	        			removed = removed + 1
	            		# cat("Removing edge ",colnames(adj.matrix)[curr.edge.i]," to ",colnames(adj.matrix)[curr.edge.j],"\n");
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
    		
    		# expanded matrix to be considered in removing the loops
    		expansion = hypotheses.expansion(input_matrix=adj.matrix,map=hypotheses$hstructure,hidden_and=F,expand=T,skip.disconnected=F);
    		
        # load the igraph library required for the loop detection
    	if (!require(igraph)) {
    		install.packages('igraph', dependencies = TRUE);
      		library(igraph);
    	}
    	
        for(i in 1:length(ordered.edges)) {
        	
            # consider the edge i-->j
            curr.edge = ordered.edges[[i]];
            curr.edge.i = curr.edge[1,1];
            curr.edge.j = curr.edge[2,1];
            
            # resolve the mapping from the adj.matrix to the expanded one both for curr.edge.i and curr.edge.j
            if(colnames(adj.matrix)[curr.edge.i]%in%expansion[[2]]) {
            		curr.edge.i.exp = which(colnames(expansion[[1]])%in%names(expansion[[2]])[which(expansion[[2]]%in%colnames(adj.matrix)[curr.edge.i])]);
            }
            else {
            		curr.edge.i.exp = which(colnames(expansion[[1]])%in%colnames(adj.matrix)[curr.edge.i]);
            }
            if(colnames(adj.matrix)[curr.edge.j]%in%expansion[[2]]) {
            		curr.edge.j.exp = which(colnames(expansion[[1]])%in%names(expansion[[2]])[which(expansion[[2]]%in%colnames(adj.matrix)[curr.edge.j])]);
            }
            else {
            		curr.edge.j.exp = which(colnames(expansion[[1]])%in%colnames(adj.matrix)[curr.edge.j]);
            }
            
            # search for loops between curr.edge.i and curr.edge.j
            curr.graph = graph.adjacency(expansion[[1]], mode="directed");
            	is.path = length(unlist(get.shortest.paths(curr.graph, curr.edge.j.exp, curr.edge.i.exp)$vpath));
            
            # if there is a path between the two nodes, remove edge i --> j
            if(is.path>0) {
            		removed = removed + 1
            		# cat("Removing edge ",colnames(adj.matrix)[curr.edge.i]," to ",colnames(adj.matrix)[curr.edge.j],"\n");
            		expansion[[1]][curr.edge.i.exp,curr.edge.j.exp] = 0;
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
