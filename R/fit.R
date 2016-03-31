# perform the likelihood fit
lregfit <- function(data,
                    adj.matrix,
                    adj.matrix.fit,
                    regularization,
                    command){

    if (regularization == "no_reg") {
        return(adj.matrix)
    }

    # create the blacklist based on the prima facie topology and the tree-structure assumption
    cont = 0
    parent = -1
    child = -1
    for (i in 1:nrow(adj.matrix)) {
        for (j in 1:ncol(adj.matrix)) {
            if(i != j) {
                if (adj.matrix[i,j] == 0) {
                    # [i,j] refers to causation i --> j
                    cont = cont + 1
                    if (cont == 1) {
                        parent = toString(i)
                        child = toString(j)
                    } else {
                        parent = c(parent, toString(i))
                        child = c(child, toString(j))
                    }
                }
            }
        }
    }

    # perform the reconstruction by likelihood fit with regularization
    # either the hill climbing or the tabu search is used as the mathematical optimization technique
    if (cont > 0) {
        blacklist = data.frame(from = parent,to = child)
        if (command == "hc") {
            my.net = hc(data,score = regularization, blacklist = blacklist)
        } else if (command == "tabu") {
            my.net = tabu(data,score = regularization, blacklist = blacklist)
        }
    } else {
        if (command == "hc") {
            my.net = hc(data, score = regularization)
        } else if (command == "tabu") {
            my.net = tabu(data, score = regularization)
        }
    }
    my.arcs = my.net$arcs
    
    # build the adjacency matrix of the reconstructed topology
    if (length(nrow(my.arcs)) > 0 && nrow(my.arcs) > 0) {
        for (i in 1:nrow(my.arcs)) {
            # [i,j] refers to causation i --> j
            adj.matrix.fit[as.numeric(my.arcs[i,1]), 
                           as.numeric(my.arcs[i,2])] = 1
        }
    }

    return(adj.matrix.fit)
}