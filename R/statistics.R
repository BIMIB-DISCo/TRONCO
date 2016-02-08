#### TRONCO: a tool for TRanslational ONCOlogy
####
#### Copyright (c) 2015-2016, Marco Antoniotti, Giulio Caravagna, Luca De Sano,
#### Alex Graudenzi, Ilya Korsunsky, Mattia Longoni, Loes Olde Loohuis,
#### Giancarlo Mauri, Bud Mishra and Daniele Ramazzotti.
####
#### All rights reserved. This program and the accompanying materials
#### are made available under the terms of the GNU GPL v3.0
#### which accompanies this distribution.


#' Convert a TRONCO object in a Bnlearn network.
#' @title as.bnlearn.network
#'
#' @examples
#' data(test_model)
#' as.bnlearn.network(test_model)
#'
#' @param data A reconstructed model (the output of tronco.capri or tronco.caprese)
#' @param regularization The name of the selected regularization (default: "bic")
#' @param makeValid Transform the bootstrapped data into a valid 2-categories input data
#' @export as.bnlearn.network
#'
as.bnlearn.network <- function(obj, regularization = "bic", makeValid = TRUE) {

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

    if (makeValid) {
        for (i in 1:ncol(genotypes)) {
            if (sum(genotypes[, i]) == 0) {
                genotypes[sample(1:nrow(genotypes), size=1), i] = 1;
            } else if (sum(genotypes[, i]) == nrow(genotypes)) {
                genotypes[sample(1:nrow(genotypes), size=1), i] = 0;
            }
        }
    }

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


#' Perform a k-fold cross-validation using the function bn.cv
#' to estimate the entropy loss.
#' @title tronco.kfold.eloss
#'
#' @examples
#' data(test_model)
#' tronco.kfold.eloss(test_model)
#'
#' @param data A reconstructed model (the output of tronco.capri or tronco.caprese)
#' @param regularization The name of the selected regularization (default: "bic")
#' @param runs a positive integer number, the number of times cross-validation will be run
#' @param k a positive integer number, the number of groups into which the data will be split
#' @importFrom bnlearn bn.cv
#' @export tronco.kfold.eloss
#'
tronco.kfold.eloss = function(x, 
                              regularization = as.parameters(x)$regularization,
                              runs = 10,
                              k = 10) {   

    ## Check if there is a reconstructed model.

    if (!has.model(x)) {
        stop('The input TRONCO object does not contain a model, you should first do that -- won\'t perform cross-validation!')
    }

    ## Check if the reconstruction has been made with CAPRI

    if (x$parameters$algorithm != 'CAPRI') {
        stop('The model contained in the input TRONCO object has not been reconstructed with CAPRI,  -- won\'t perform cross-validation!')
    }

    ## Check if the selected regularization is used in the model.

    if (!"kfold" %in% names(x)) {
        x$kfold = NULL
    }

    for (reg in regularization) {
        if (!reg %in% as.parameters(x)$regularization) {
            stop(paste(reg, " was not used to infer the input TRONCO object -- won\'t perform cross-validation!"))
        }

        ## Get bnlearn network.
        bn = as.bnlearn.network(x, reg)
        bndata = bn$data
        bnnet = bn$net

        ## Calculating the eloss with bn.cv
        cat('Calculating entropy loss with k-fold cross-validation [ k =', k,
            '| runs =', runs, '| regularizer =', reg, '] ... ')
        bn.kcv.list = bn.cv(bndata, bnnet, loss = 'logl', runs = runs, k = k)

        losses = NULL
        for(i in 1:length(bn.kcv.list)) {
            losses = c(losses, attr(bn.kcv.list[[i]], "mean"))
        }
        if (any(is.na(losses))) {
            warning("Some folds returned NA")
        }
        eloss = losses[!is.na(losses)]
    
        meanll = mean(eloss)
        ll = get(reg, x$model)$logLik
        ratio = meanll / abs(ll) * 100
    
        cat(' DONE\n')
        cat('  Model logLik =', ll, '\n')
        cat('  Mean   eloss =', meanll,' | ', ratio,'% \n')
        cat('  Stdev  eloss = ', sd(eloss) ,'\n')
    
        x$kfold[[reg]]$bn.kcv.list = bn.kcv.list
        x$kfold[[reg]]$eloss = eloss
    }
    return(x)
}


#' Perform a k-fold cross-validation using the function bn.cv
#' and scan every node to estimate its prediction error. 
#' @title tronco.kfold.prederr
#'
#' @examples
#' data(test_model)
#' tronco.kfold.prederr(test_model)
#'
#' @param x A reconstructed model (the output of tronco.capri)
#' @param regularization The name of the selected regularization (default: "bic")
#' @param events a list of event 
#' @param runs a positive integer number, the number of times cross-validation will be run
#' @param k a positive integer number, the number of groups into which the data will be split
#' @param cores.ratio Percentage of cores to use. coresRate * (numCores - 1)
#' @importFrom bnlearn bn.cv
#' @export tronco.kfold.prederr
#'
tronco.kfold.prederr <- function(x,
                                 regularization = as.parameters(x)$regularization,
                                 events = as.events(x),
                                 runs = 10,
                                 k = 10,
                                 cores.ratio = 1,
                                 verbose = FALSE) {

    ## Check if there is a reconstructed model.

    if(!has.model(x)) {
        stop('This object does not have a model.')
    }

    ## Check if the reconstruction has been made with CAPRI

    if (x$parameters$algorithm != 'CAPRI') {
        stop('The model contained in the input TRONCO object has not been reconstructed with CAPRI,  -- won\'t perform cross-validation!')
    }

    ## Integrity check over nodes.
    adj.matrix = as.adj.matrix(x,
                               events = events,
                               models = regularization,
                               type = 'fit')
  
    events = apply(events, 1, function(z){paste(z, collapse = ' ')})

    if (!"kfold" %in% names(x)) {
        x$kfold = NULL
    }

    cores = as.integer(cores.ratio * (detectCores() - 1))
    if (cores < 1) {
        cores = 1
    }
    if (verbose) {
        cl = makeCluster(cores, outfile = '')
    } else {
        cl = makeCluster(cores)
    }
    registerDoParallel(cl)
    cat('*** Using', cores, 'cores via "parallel" \n')


    for (reg in regularization) {
        
        ## Check if the selected regularization is used in the model.
        
        if (!reg %in% as.parameters(x)$regularization) {
            stop(paste(reg, " was not used to infer the input TRONCO object -- won\'t perform cross-validation!"))
        }

        ## Get bnlearn network and the adj matrix.

        bn = as.bnlearn.network(x, reg)
        bndata = bn$data
        bnnet = bn$net

        adj.matrix = get(reg, as.adj.matrix(x))
        adj.matrix = keysToNames(x, adj.matrix)
        names(colnames(adj.matrix)) = NULL
        names(rownames(adj.matrix)) = NULL
  
        pred = list()

        ## Perform the estimation of the prediction error. 
        
        cat('Scanning', length(events), 'nodes for their prediction error (all parents). Regularizer: ', reg, '\n')

        

        r = foreach(i = 1:length(events), .inorder = TRUE) %dopar% {
            cat('\tprocessing ', events[i], '\n')

            comp = bn.cv(bndata,
                         bnnet, 
                         loss = 'pred', 
                         loss.args = list(target = events[i]),
                         runs = runs,
                         k = k)
            
            res = NULL
            for(i in 1:runs) {
                res = c(res, attributes(comp[[i]])$mean)
            }
 
            res  
        }
        cat("*** Reducing results\n")
        names(r) = events      
        x$kfold[[reg]]$prederr = r
    }

    stopCluster(cl)

    return(x)
}

#' Perform a k-fold cross-validation using the function bn.cv
#' and scan every node to estimate its posterior classification error. 
#' @title tronco.kfold.posterr
#'
#' @examples
#' data(test_model)
#' tronco.kfold.posterr(test_model)
#'
#' @param x A reconstructed model (the output of tronco.capri)
#' @param regularization The name of the selected regularization (default: "bic")
#' @param events a list of event 
#' @param runs a positive integer number, the number of times cross-validation will be run
#' @param k a positive integer number, the number of groups into which the data will be split
#' @importFrom bnlearn bn.cv
#' @export tronco.kfold.posterr
#'
tronco.kfold.posterr <- function(x,
                                 regularization = as.parameters(x)$regularization,
                                 events = as.events(x),
                                 runs = 10,
                                 k = 10) {

    ## Check if there is a reconstructed model.

    if(!has.model(x)) {
        stop('This object does not have a model.')
    }

    ## Check if the reconstruction has been made with CAPRI

    if (x$parameters$algorithm != 'CAPRI') {
        stop('The model contained in the input TRONCO object has not been reconstructed with CAPRI,  -- won\'t perform cross-validation!')
    }

    ## Integrity check over nodes.
    adj.matrix = as.adj.matrix(x,
                               events = events,
                               models = regularization,
                               type = 'fit')
  
    events = apply(events, 1, function(z){paste(z, collapse = ' ')})

    if (!"kfold" %in% names(x)) {
        x$kfold = NULL
    }


    for (reg in regularization) {
        
        ## Check if the selected regularization is used in the model.
        
        if (!reg %in% as.parameters(x)$regularization) {
            stop(paste(reg, " was not used to infer the input TRONCO object -- won\'t perform cross-validation!"))
        }

        ## Get bnlearn network and the adj matrix.

        bn = as.bnlearn.network(x, reg)
        bndata = bn$data
        bnnet = bn$net

        adj.matrix = get(reg, as.adj.matrix(x))
        adj.matrix = keysToNames(x, adj.matrix)
        names(colnames(adj.matrix)) = NULL
        names(rownames(adj.matrix)) = NULL
        posterr.adj.matrix = array(list(NA), c(nrow(adj.matrix), ncol(adj.matrix)))
        colnames(posterr.adj.matrix) = colnames(adj.matrix)
        rownames(posterr.adj.matrix) = rownames(adj.matrix)
  
        pred = list()

        ## Perform the estimation of the prediction error. 
        
        cat('Scanning', sum(adj.matrix == 1), 'edges for posterior classification error. Regularizer:', reg, '\n')
        for (event in events) {
            if (any(adj.matrix[,event,drop=F])) {
                cat('\nTarget: ', event)
            }
            for (pre in rownames(posterr.adj.matrix)) {
                if (adj.matrix[pre,event] == 1) {

                    cat('\n\t from: ', pre, ': ')

                    comp = bn.cv(bndata,
                                 bnnet, 
                                 loss = 'pred-lw', 
                                 loss.args = list(target = event, from = pre),
                                 runs = runs,
                                 k = k)
                    
                    res = NULL
                    for(i in 1:runs) {
                        res = c(res, attributes(comp[[i]])$mean)
                    }
                    cat(mean(res), ' (', sd(res), ')')
                    posterr.adj.matrix[[pre,event]] = res  
                }
            }
        }
        x$kfold[[reg]]$posterr = posterr.adj.matrix
        cat('\n')
    }
    return(x)
}
