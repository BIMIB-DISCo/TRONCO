#### TRONCO: a tool for TRanslational ONCOlogy
####
#### Copyright (c) 2015-2016, Marco Antoniotti, Giulio Caravagna, Luca De Sano,
#### Alex Graudenzi, Ilya Korsunsky, Mattia Longoni, Loes Olde Loohuis,
#### Giancarlo Mauri, Bud Mishra and Daniele Ramazzotti.
####
#### All rights reserved. This program and the accompanying materials
#### are made available under the terms of the GNU GPL v3.0
#### which accompanies this distribution.


# Convert a TRONCO object in a Bnlearn network.
# @title as.bnlearn.network
#
# @examples
# data(test_model)
# as.bnlearn.network(test_model)
#
# @param x A reconstructed model (the output of tronco.capri or tronco.caprese)
# @param model The name of the selected regularization
# @param makeValid Transform the bootstrapped data into a valid 2-categories input data
# @export as.bnlearn.network
# @importFrom bnlearn empty.graph set.arc
#
as.bnlearn.network <- function(x, 
                               model = names(as.models(x))[1], 
                               make.valid = TRUE) {

    ## Check if there is a reconstructed model.

    if(!has.model(x)) {
        stop('Input doesn\'t have a TRONCO object inside.')
    }

    ## Check if the selected regularization is used in the model.

    if (!model %in% names(as.models(x))) {
        stop(paste(model, " was not used to build the input TRONCO model!"))
    }

    ## Get genotypes and data.

    genotypes = as.genotypes(x)
    genotypes = as.matrix(genotypes)
    genotypes = keysToNames(x, genotypes)
    names(colnames(genotypes)) = NULL

    df = as.categorical.dataset(genotypes, make.valid = make.valid)

    adj.matrix = get(model, as.adj.matrix(x, models = model))
    adj.matrix = keysToNames(x, adj.matrix)
    names(colnames(adj.matrix)) = NULL
    names(rownames(adj.matrix)) = NULL
    
    bayes.net = NULL
    bayes.net$data = df
        
    ## Create the Bayesian Network of the fitted model.
    bayes.net$net = empty.graph(colnames(genotypes))
    for (i in 1:nrow(adj.matrix)) {
        for(j in 1:ncol(adj.matrix)) {
            if(adj.matrix[i,j] == 1) {
                bayes.net$net = set.arc(
                    bayes.net$net, 
                    from = colnames(genotypes)[i], 
                    to = colnames(genotypes)[j])
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
#' tronco.kfold.eloss(test_model, k = 2, runs = 2)
#'
#' @param x A reconstructed model (the output of tronco.capri or tronco.caprese)
#' @param models The names of the selected regularizers (bic, aic or caprese)
#' @param runs a positive integer number, the number of times cross-validation will be run
#' @param k a positive integer number, the number of groups into which the data will be split
#' @param silent A parameter to disable/enable verbose messages.
#' @importFrom bnlearn bn.cv empty.graph set.arc
#' @importFrom stats sd
#' @export tronco.kfold.eloss
#'
tronco.kfold.eloss = function(x, 
                              models = names(as.models(x)),
                              runs = 10,
                              k = 10, 
                              silent = FALSE) {   

    ## Check if there is a reconstructed model.

    if (!has.model(x)) {
        stop('The input TRONCO object does not contain a model, you should first do that -- won\'t perform cross-validation!')
    }

    ## Check if the selected regularization is used in the model.

    if (!"kfold" %in% names(x)) {
        x$kfold = NULL
    }

    for (model in models) {
        if (!model %in% names(as.models(x))) {
            stop(paste(model, " was not used to infer the input TRONCO object -- won\'t perform cross-validation!"))
        }

        ## Get bnlearn network.
        bn = as.bnlearn.network(x, model)
        bndata = bn$data
        bnnet = bn$net

        ## Calculating the eloss with bn.cv
        if (!silent) {
            cat('Calculating entropy loss with k-fold cross-validation \n\t[ k =', k,
                '| runs =', runs, '| regularizer =', model, '] ... ')
        }

        ## Scutari fix
        
        bn.kcv.list = bn.cv(bndata, 
                            bnnet,
                            loss = 'logl',
                            runs = runs,
                            k = k,
                            fit = "bayes",
                            fit.args = list(iss = 1))

        losses = NULL

        for(i in 1:length(bn.kcv.list)) {
            losses = c(losses, attr(bn.kcv.list[[i]], "mean"))
        }
        if (any(is.na(losses))) {
            warning("Some folds returned NA")
        }
        eloss = losses[!is.na(losses)]
    
        meanll = mean(eloss)
        ll = get(model, x$model)$logLik
        ratio = meanll / abs(ll) * 100
    
        if (!silent) {
            cat(' DONE\n')
            cat('  Model logLik =', ll, '\n')
            cat('  Mean   eloss =', meanll,' | ', ratio,'% \n')
            cat('  Stdev  eloss = ', sd(eloss) ,'\n')
        }
    
        x$kfold[[model]]$bn.kcv.list = bn.kcv.list
        x$kfold[[model]]$eloss = eloss
    }
    return(x)
}


#' Perform a k-fold cross-validation using the function bn.cv
#' and scan every node to estimate its prediction error. 
#' @title tronco.kfold.prederr
#'
#' @examples
#' data(test_model)
#' tronco.kfold.prederr(test_model, k = 2, runs = 2)
#'
#' @param x A reconstructed model (the output of tronco.capri)
#' @param models The names of the selected regularizers (bic, aic or caprese)
#' @param events a list of event 
#' @param runs a positive integer number, the number of times cross-validation will be run
#' @param k a positive integer number, the number of groups into which the data will be split
#' @param cores.ratio Percentage of cores to use. coresRate * (numCores - 1)
#' @param silent A parameter to disable/enable verbose messages.
#' @importFrom bnlearn bn.cv empty.graph set.arc
#' @importFrom doParallel registerDoParallel  
#' @importFrom foreach foreach %dopar%
#' @importFrom iterators icount
#' @importFrom parallel stopCluster makeCluster detectCores
#' @importFrom stats sd
#' @export tronco.kfold.prederr
#'
tronco.kfold.prederr <- function(x,
                                 models = names(as.models(x)),
                                 events = as.events(x),
                                 runs = 10,
                                 k = 10,
                                 cores.ratio = 1,
                                 silent = FALSE) {

    ## Check if there is a reconstructed model.

    if(!has.model(x)) {
        stop('This object does not have a model.')
    }

    ## Integrity check over nodes.
    as.adj.matrix(x,
                  events = events,
                  models = models,
                  type = 'fit')
  
    events = apply(events, 1, function(z){paste(z, collapse = ' ')})

    if (!"kfold" %in% names(x)) {
        x$kfold = NULL
    }

    cores = as.integer(cores.ratio * (detectCores() - 1))
    if (cores < 1) {
        cores = 1
    }

    cl = makeCluster(cores)
    registerDoParallel(cl)

    if (!silent) {
        cat('*** Using', cores, 'cores via "parallel" \n')
    }


    for (model in models) {
        
        ## Check if the selected regularization is used in the model.
        
        if (!model %in% names(as.models(x))) {
            stop(paste(model, " was not used to infer the input TRONCO object -- won\'t perform cross-validation!"))
        }

        ## Get bnlearn network and the adj matrix.

        bn = as.bnlearn.network(x, model)
        bndata = bn$data
        bnnet = bn$net

        adj.matrix = get(model, as.adj.matrix(x))
        adj.matrix = keysToNames(x, adj.matrix)
        names(colnames(adj.matrix)) = NULL
        names(rownames(adj.matrix)) = NULL
  
        pred = list()

        ## Perform the estimation of the prediction error. 
        if (!silent) {
            cat('\tScanning',
                length(events),
                'nodes for their prediction error (all parents). \n\tRegularizer: ',
                model,
                '\n')
        }
        

        r = foreach(i = 1:length(events), .inorder = TRUE) %dopar% {

            ## Scutari fix

            event = events[i]

            comp = bn.cv(bndata,
                         bnnet, 
                         loss = 'pred', 
                         loss.args = list(target = event),
                         runs = runs,
                         k = k,
                         fit = "bayes",
                         fit.args = list(iss = 1))
            
            res = NULL
            for(i in 1:runs) {
                res = c(res, attributes(comp[[i]])$mean)
            }            
            res  
        }
        if (!silent) {
            cat("*** Reducing results\n")
        }
        names(r) = events      
        x$kfold[[model]]$prederr = r
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
#' tronco.kfold.posterr(test_model, k = 2, runs = 2)
#'
#' @param x A reconstructed model (the output of tronco.capri)
#' @param models The names of the selected regularizers (bic, aic or caprese)
#' @param events a list of event 
#' @param runs a positive integer number, the number of times cross-validation will be run
#' @param k a positive integer number, the number of groups into which the data will be split
#' @param cores.ratio Percentage of cores to use. coresRate * (numCores - 1)
#' @param silent A parameter to disable/enable verbose messages.
#' @importFrom bnlearn bn.cv empty.graph set.arc
#' @importFrom stats sd
#' @export tronco.kfold.posterr
#'
tronco.kfold.posterr <- function(x,
                                 models = names(as.models(x)),
                                 events = as.events(x),
                                 runs = 10,
                                 k = 10,
                                 cores.ratio = 1,
                                 silent = FALSE) {

    ## Check if there is a reconstructed model.

    if(!has.model(x)) {
        stop('This object does not have a model.')
    }

    ## Integrity check over nodes.
    as.adj.matrix(x,
                  events = events,
                  models = models,
                  type = 'fit')
  
    events = apply(events, 1, function(z){paste(z, collapse = ' ')})

    if (!"kfold" %in% names(x)) {
        x$kfold = NULL
    }


    cores = as.integer(cores.ratio * (detectCores() - 1))
    if (cores < 1) {
        cores = 1
    }

    cl = makeCluster(cores)

    registerDoParallel(cl)
    if (!silent) {
        cat('*** Using', cores, 'cores via "parallel" \n')
    }

    for (model in models) {
        
        ## Check if the selected regularization is used in the model.
        
        if (!model %in% names(as.models(x))) {
            stop(paste(model, " was not used to infer the input TRONCO object -- won\'t perform cross-validation!"))
        }

        ## Get bnlearn network and the adj matrix.

        bn = as.bnlearn.network(x, model)
        bndata = bn$data
        bnnet = bn$net

        adj.matrix = get(model, as.adj.matrix(x))
        adj.matrix = keysToNames(x, adj.matrix)
        names(colnames(adj.matrix)) = NULL
        names(rownames(adj.matrix)) = NULL
        


        ## Perform the estimation of the prediction error. 
        
        if (!silent) {
            cat('\tScanning',
                sum(adj.matrix == 1),
                'edges for posterior classification error. \n\tRegularizer:',
                model,
                '\n')
        }    

        r = foreach(i = 1:length(events), .inorder = TRUE, .combine = cbind) %dopar% {
            event = events[i]
            posterr.adj.col = array(list(NA), c(nrow(adj.matrix), 1))
            colnames(posterr.adj.col) = event
            rownames(posterr.adj.col) = rownames(adj.matrix)
                        
            for (pre in rownames(posterr.adj.col)) {
                if (adj.matrix[pre,event] == 1) {

                    ## Scutari fix

                    comp = bn.cv(bndata,
                                 bnnet, 
                                 loss = 'pred-lw', 
                                 loss.args = list(target = event, from = pre),
                                 runs = runs,
                                 k = k,
                                 fit = "bayes",
                                 fit.args = list(iss = 1))
                    
                    res = NULL
                    for(i in 1:runs) {
                        res = c(res, attributes(comp[[i]])$mean)
                    }
                    posterr.adj.col[[pre,event]] = res  
                }
            }
            posterr.adj.col
        }
        if (!silent) {
            cat("*** Reducing results\n")
        }
        x$kfold[[model]]$posterr = r
    }

    stopCluster(cl)
    return(x)
}
