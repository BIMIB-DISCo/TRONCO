#### TRONCO: a tool for TRanslational ONCOlogy
####
#### Copyright (c) 2015-2016, Marco Antoniotti, Giulio Caravagna, Luca De Sano,
#### Alex Graudenzi, Giancarlo Mauri, Bud Mishra and Daniele Ramazzotti.
####
#### All rights reserved. This program and the accompanying materials
#### are made available under the terms of the GNU GPL v3.0
#### which accompanies this distribution.


# perform the likelihood fit
lregfit <- function(data,
                    adj.matrix,
                    adj.matrix.fit,
                    regularization,
                    command){

    if (regularization == "no_reg") {
        return(adj.matrix)
    }


    ## Create the blacklist based on the prima facie topology
    ## and the tree-structure assumption
    cont = 0
    parent = -1
    child = -1

    for (i in rownames(adj.matrix)) {
        for (j in colnames(adj.matrix)) {
            if(i != j) {
                if (adj.matrix[i,j] == 0) {
                    # [i,j] refers to causation i --> j
                    cont = cont + 1
                    if (cont == 1) {
                        parent = i
                        child = j
                    } else {
                        parent = c(parent, i)
                        child = c(child, j)
                    }
                }
            }
        }
    }

    ## perform the reconstruction by likelihood fit with regularization
    ## either the hill climbing or the tabu search is used as the 
    ## mathematical optimization technique
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
    
    ## build the adjacency matrix of the reconstructed topology
    if (length(nrow(my.arcs)) > 0 && nrow(my.arcs) > 0) {
        for (i in 1:nrow(my.arcs)) {
            # [i,j] refers to causation i --> j
            adj.matrix.fit[my.arcs[i,1], my.arcs[i,2]] = 1
        }
    }

    return(adj.matrix.fit)
}

#### end of file -- fit.R
