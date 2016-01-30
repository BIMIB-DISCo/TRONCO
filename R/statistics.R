#### TRONCO: a tool for TRanslational ONCOlogy
####
#### Copyright (c) 2015-2016, Marco Antoniotti, Giulio Caravagna, Luca De Sano,
#### Alex Graudenzi, Ilya Korsunsky, Mattia Longoni, Loes Olde Loohuis,
#### Giancarlo Mauri, Bud Mishra and Daniele Ramazzotti.
####
#### All rights reserved. This program and the accompanying materials
#### are made available under the terms of the GNU GPL v3.0
#### which accompanies this distribution.


#' Convert a TRONCO object in a Bnlean network
#' @title as.bnlearn.network
#'
#' @examples
#' data(test_model)
#' as.bnlearn.network(test_model)
#'
#' @param data A reconstructed model (the output of tronco.capri or tronco.caprese)
#' @param regularization The name of the selected regularization (default: "bic")
#' @export as.bnlearn.network
#'
as.bnlearn.network <- function(obj, regularization = "bic") {

    ## Check if there is a reconstructed model.

    if(!has.model(obj)) {
        stop('Input obj doesn\'t have a TRONCO object inside.')
    }

    ## Check if the selected regularization is used in the model.

    if (!regularization %in% names(obj$model)) {
        stop(paste(regularization, " was not used to build the input TRONCO model!"))
    }

    ## Get genotypes and data.

    genotypes = as.genotypes(obj)
    genotypes = keysToNames(obj, genotypes)
    names(colnames(genotypes)) = NULL

    adj.matrix = get(regularization, as.adj.matrix(obj))
    adj.matrix = keysToNames(obj, adj.matrix)
    names(colnames(adj.matrix)) = NULL
    names(rownames(adj.matrix)) = NULL
    
    bayes.net = NULL
            
    ## Create a categorical data frame from the dataset.
    df = array("missing",c(nrow(genotypes),ncol(genotypes)))
    for (i in 1:nrow(genotypes)) {
        for (j in 1:ncol(genotypes)) {
            if(genotypes[i,j]==1) {
                df[i,j] = "observed"
            }
        }
    }
    df = as.data.frame(df)
    my.names = names(df)
    
    for (i in 1:length(my.names)) {
        my.names[i] = toString(i)
    }
    colnames(df) = colnames(genotypes)
    bayes.net$data = df
        
    ## Create the Bayesian Network of the fitted model.
    bayes.net$net = empty.graph(colnames(genotypes))
    for (i in 1:nrow(adj.matrix)) {
        for(j in 1:ncol(adj.matrix)) {
            if(adj.matrix[i,j]==1) {
                bayes.net$net = set.arc(
                    bayes.net$net, 
                    from=colnames(genotypes)[i], 
                    to=colnames(genotypes)[j])
            }
        }
    }
    
        
    return(bayes.net) 
}


#' Perform a k-fold cross-validation using the function bn.cv.
#' @title stat.eloss
#'
#' @examples
#' data(test_model)
#' stat.eloss(test_model)
#'
#' @param data A reconstructed model (the output of tronco.capri or tronco.caprese)
#' @param regularization The name of the selected regularization (default: "bic")
#' @param runs a positive integer number, the number of times cross-validation will be run
#' @param k a positive integer number, the number of groups into which the data will be split
#' @importFrom bnlearn bn.cv
#' @export stat.eloss
#'
stat.eloss = function(data, regularization = "bic", runs = 10, k = 10) {   
    
    ## Check if there is a reconstructed model.

    if (!has.model(data)) {
        stop('This dataset doesn\'t have.')
    }

    ## Check if the selected regularization is used in the model.

    if (!regularization %in% names(data$model)) {
        stop(paste(regularization, "not in model"))
    }

    ## Get bnlearn network.

    bn = as.bnlearn.network(data, regularization)
    bndata = bn$data
    bnnet = bn$net

    ## Calculating the eloss with bn.cv

    eloss = NULL
    cat('Entropy loss ...')
    bn.kcv.list = bn.cv(bndata, bnnet, loss = 'logl', runs = runs, k = k)
    eloss$bn.kcv.list = bn.kcv.list

    losses = NULL
    for(i in 1:length(bn.kcv.list)) {
        losses = c(losses, attr(bn.kcv.list[[i]], "mean"))
    }
    if (any(is.na(losses))) {
        warning("Some run returned NA")
    }
    eloss$value = losses[!is.na(losses)]
    message(' DONE')
    return(eloss)
}


#' Perform a k-fold cross-validation (with k = 10) using the function bn.cv
#' and scan every node to estimate its prediction error. 
#' @title stat.prederr
#'
#' @examples
#' data(test_model)
#' stat.prederr(test_model)
#'
#' @param data A reconstructed model (the output of tronco.capri or tronco.caprese)
#' @param regularization The name of the selected regularization (default: "bic")
#' @param nodes a list of event 
#' @param runs a positive integer number, the number of times cross-validation will be run
#' @importFrom bnlearn bn.cv
#' @export stat.prederr
#'
stat.prederr <- function(data,
                         regularization = "bic", 
                         nodes = as.events(data, keysToNames = TRUE),
                         runs = 10) {

    ## Check if there is a reconstructed model.

    if(!has.model(data)) {
        stop('This dataset doesn\'t have.')
    }

    ## Check if the selected regularization is used in the model.

    if (!regularization %in% names(data$model)) {
        stop(paste(regularization, "not in model"))
    }

    ## Get bnlearn network and the adj matrix.

    bn = as.bnlearn.network(data, regularization)
    bndata = bn$data
    bnnet = bn$net

    adj.matrix = get(regularization, as.adj.matrix(data))
    adj.matrix = keysToNames(data, adj.matrix)
    names(colnames(adj.matrix)) = NULL
    names(rownames(adj.matrix)) = NULL

    ## Integrity check over nodes.
    for(node in nodes) {
        if(!node %in% rownames(adj.matrix)) {
            stop(paste("Invalid node found: ", node))
        }
    }

    pred = list()

    ## Perform the estimation of the prediction error. 
    
    cat('Scanning', length(nodes), 'nodes for their prediction error.\n')
    for (i in 1:length(nodes)) {
        cat(i, 'Prediction error (parents):', nodes[i], '...')

        comp = bn.cv(bndata,
                     bnnet, 
                     loss = 'pred', 
                     loss.args = list(target = nodes[i]),
                     runs = runs)
        
        res = NULL
        for(i in 1:runs) {
            res = c(res, attributes(comp[[i]])$mean)
        }

        pred = append(pred, list(res))
        message(' DONE')    
    }
    names(pred) = nodes
    return(pred)
}
