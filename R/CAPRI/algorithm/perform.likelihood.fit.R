#### perform.likelihood.fit.R
####
#### TRONCO: a tool for TRanslational ONCOlogy
####
#### See the files COPYING and LICENSE for copyright and licensing
#### information.


# reconstruct the best causal topology by likelihood fit
# INPUT:
# dataset: a valid dataset
# adj.matrix: the adjacency matrix of the prima facie causes
# command: type of search, either hill climbing (hc) or tabu (tabu)
# regularization: regularization term to be used in the likelihood fit
# RETURN:
# topology: the adjacency matrix of both the prima facie and causal topologies
"perform.likelihood.fit" <-
function( dataset, adj.matrix, command, regularization ) {
	
    # load the bnlearn library required for the likelihood fit with bic
    if (!require(bnlearn)) {
    		install.packages('bnlearn', dependencies = TRUE);
      	library(bnlearn);
    }
    
    # adjacency matrix of the topology reconstructed by likelihood fit
    adj.matrix.fit = array(0,c(nrow(adj.matrix),ncol(adj.matrix)));

    # create a categorical data frame from the dataset
    data = array("missing",c(nrow(dataset),ncol(dataset)));
    for (i in 1:nrow(dataset)) {
        for (j in 1:ncol(dataset)) {
            if(dataset[i,j]==1) {
                data[i,j] = "observed";
            }
        }
    }
    data = as.data.frame(data);
    my.names = names(data);
    for (i in 1:length(my.names)) {
        my.names[i] = toString(i);
    }
    names(data) = my.names;
    
    # create the blacklist based on the prima facie topology
    cont = 0;
    parent = -1;
    child = -1;
    for (i in 1:nrow(adj.matrix)) {
        for (j in 1:ncol(adj.matrix)) {
            if(i!=j) {
                if(adj.matrix[i,j]==0) {
                    # [i,j] refers to causation i --> j
                    cont = cont + 1;
                    if(cont==1) {
                        parent = toString(i);
                        child = toString(j);
                    }
                    else {
                        parent = c(parent,toString(i));
                        child = c(child,toString(j));
                    }
                }
            }
        }
    }

    # perform the reconstruction by likelihood fit with regularization
    # either the hill climbing or the tabu search is used as the mathematical optimization technique
    
    # cat('Performing likelihood-fit with regularization:', regularization, '(bnlearn)\n');
    # cat('Heuristic search method:', command, '(bnlearn)\n');
    
    if(cont>0) {
        blacklist = data.frame(from = parent,to = child);
        if(command=="hc") {
        		my.net = hc(data,score= regularization,blacklist=blacklist);
        }
        else if(command=="tabu") {
        		my.net = tabu(data,score= regularization,blacklist=blacklist);
        }
    }
    else {
    		if(command=="hc") {
        		my.net = hc(data,score= regularization);
        }
        else if(command=="tabu") {
        		my.net = tabu(data,score= regularization);
        }
    }
    my.arcs = my.net$arcs;
    
    # build the adjacency matrix of the reconstructed topology
    if(length(nrow(my.arcs))>0 && nrow(my.arcs)>0) {
        for (i in 1:nrow(my.arcs)) {
            # [i,j] refers to causation i --> j
            adj.matrix.fit[as.numeric(my.arcs[i,1]),as.numeric(my.arcs[i,2])] = 1;
        }
    }
    
    # save the results and return them
    adj.matrix = list(adj.matrix.pf=adj.matrix,adj.matrix.fit=adj.matrix.fit);
    topology = list(adj.matrix=adj.matrix);
    return(topology);

}

#### end of file -- perform.likelihood.fit.R
