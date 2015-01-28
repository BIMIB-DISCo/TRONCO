#### perform.likelihood.fit.R
####
#### TRONCO: a tool for TRanslational ONCOlogy
####
#### See the files COPYING and LICENSE for copyright and licensing
#### information.


#reconstruct the best minimal causal topology by likelihood fit
#INPUT:
#dataset: a valid dataset
#adj.matrix: the adjacency matrix of the prima facie causes
#RETURN:
#topology: the adjacency matrix of both the prima facie and causal topologies
"perform.likelihood.fit" <-
function(dataset,adj.matrix) {
    #load the bnlearn library required for the likelihood fit with bic
    require(bnlearn);
    #adjacency matrix of the topology reconstructed by likelihood fit
    adj.matrix.bic = array(0,c(nrow(adj.matrix),ncol(adj.matrix)));
    #create a categorical data frame from the dataset
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
    #create the blacklist based on the prima facie topology
    cont = 0;
    parent = -1;
    child = -1;
    for (i in 1:nrow(adj.matrix)) {
        for (j in 1:ncol(adj.matrix)) {
            if(i!=j) {
                if(adj.matrix[i,j]==0) {
                    #[i,j] refers to causation i --> j
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
    #perform the reconstruction by likelihood fit with bic
    #the hill climbing is used as the mathematical optimization technique
    if(cont>0) {
        blacklist = data.frame(from = parent,to = child);
        hc.arcs = hc(data,score="bic",blacklist=blacklist)$arcs;
    }
    else {
        hc.arcs = hc(data,score="bic")$arcs;
    }
    #make the adjacency matrix of the reconstructed topology
    if(length(nrow(hc.arcs))>0 && nrow(hc.arcs)>0) {
        for (i in 1:nrow(hc.arcs)) {
            #[i,j] refers to causation i --> j
            adj.matrix.bic[as.numeric(hc.arcs[i,1]),as.numeric(hc.arcs[i,2])] = 1;
        }
    }
    topology = list(adj.matrix.pf=adj.matrix,adj.matrix.bic=adj.matrix.bic);
    return(topology);
}

#### end of file -- perform.likelihood.fit.R
