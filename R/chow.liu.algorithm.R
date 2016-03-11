#### TRONCO: a tool for TRanslational ONCOlogy
####
#### Copyright (c) 2015-2016, Marco Antoniotti, Giulio Caravagna, Luca De Sano,
#### Alex Graudenzi, Ilya Korsunsky, Mattia Longoni, Loes Olde Loohuis,
#### Giancarlo Mauri, Bud Mishra and Daniele Ramazzotti.
####
#### All rights reserved. This program and the accompanying materials
#### are made available under the terms of the GNU GPL v3.0
#### which accompanies this distribution.


# reconstruct the best topology based on probabilistic causation and Chow Liu algorithm
# @title chow.liu.fit
# @param dataset a dataset describing a progressive phenomenon
# @param regularization regularizators to be used for the likelihood fit
# @param do.boot should I perform bootstrap? Yes if TRUE, no otherwise
# @param nboot integer number (greater than 0) of bootstrap sampling to be performed
# @param pvalue pvalue for the tests (value between 0 and 1)
# @param min.boot minimum number of bootstrapping to be performed
# @param min.stat should I keep bootstrapping untill I have nboot valid values?
# @param boot.seed seed to be used for the sampling
# @param do.estimation should I perform the estimation of the error rates and probabilities?
# @param silent should I be verbose?
# @return topology: the reconstructed tree topology
#
chow.liu.fit <- function(dataset,
                      regularization = c("bic","aic"),
                      do.boot = TRUE,
                      nboot = 100,
                      pvalue = 0.05,
                      min.boot = 3,
                      min.stat = TRUE,
                      boot.seed = NULL,
                      do.estimation = FALSE,
                      silent = FALSE ) {

    ## Start the clock to measure the execution time.
    
    ptm = proc.time();

    ## Structure with the set of valid edges
    ## I start from the complete graph, i.e., I have no prior and all
    ## the connections are possibly causal.
    
    adj.matrix = array(1, c(ncol(dataset), ncol(dataset)));
    colnames(adj.matrix) = colnames(dataset);
    rownames(adj.matrix) = colnames(dataset);

    ## The diagonal of the adjacency matrix should not be considered,
    ## i.e., no self cause is allowed.
    
    diag(adj.matrix) = 0;

    ## Check if the dataset is valid.
    
    valid.dataset = check.dataset(dataset, adj.matrix, FALSE);
    adj.matrix = valid.dataset$adj.matrix;
    invalid.events = valid.dataset$invalid.events;

    ## Reconstruct the prima facie topology
    ## Should I perform bootstrap? Yes if TRUE, no otherwise.
    
    if (do.boot == TRUE) {
        if (!silent)
            cat('*** Bootstraping selective advantage scores (prima facie).\n')
        prima.facie.parents =
            get.prima.facie.parents.do.boot(dataset,
                                            NA,
                                            nboot,
                                            pvalue,
                                            adj.matrix,
                                            min.boot,
                                            min.stat,
                                            boot.seed,
                                            silent);
    } else {
        if (!silent)
            cat('*** Computing selective advantage scores (prima facie).\n')
        prima.facie.parents =
            get.prima.facie.parents.no.boot(dataset,
                                            NA,
                                            adj.matrix,
                                            silent);
    }

    ## Add back in any connection invalid for the probability raising
    ## theory.
    
    if (length(invalid.events) > 0) {
        for (i in 1:nrow(invalid.events)) {
            prima.facie.parents$adj.matrix$adj.matrix.acyclic[invalid.events[i, "cause"],invalid.events[i, "effect"]] = 1;
            prima.facie.parents$adj.matrix$adj.matrix.cyclic[invalid.events[i, "cause"],invalid.events[i, "effect"]] = 1;
        }
    }
    adj.matrix.prima.facie =
        prima.facie.parents$adj.matrix$adj.matrix.cyclic

    ## Perform the likelihood fit with the required strategy.
    
    model = list();
    for (reg in regularization) {

        ## Perform the likelihood fit with the chosen regularization
        ## score on the prima facie topology.
        
        if (!silent)
            cat(paste0('*** Performing likelihood-fit with regularization ',reg,'.\n'))
        best.parents =
            perform.likelihood.fit.chow.liu(dataset,
                                   prima.facie.parents$adj.matrix$adj.matrix.acyclic,
                                   command,
                                   regularization = reg);

        ## Set the structure to save the conditional probabilities of
        ## the reconstructed topology.
        
        parents.pos.fit = array(list(),c(ncol(dataset),1));
        conditional.probs.fit = array(list(),c(ncol(dataset),1));

        ## Compute the conditional probabilities.
        
        for (i in 1:ncol(dataset)) {
            for (j in 1:ncol(dataset)) {
                if (i!=j && best.parents$adj.matrix$adj.matrix.fit[i, j] == 1) {
                    parents.pos.fit[j,1] =
                        list(c(unlist(parents.pos.fit[j, 1]), i));
                    conditional.probs.fit[j,1] =
                        list(c(unlist(conditional.probs.fit[j, 1]),
                               prima.facie.parents$joint.probs[i, j] /
                                   prima.facie.parents$marginal.probs[i]));
                }
            }
        }
        parents.pos.fit[unlist(lapply(parents.pos.fit, is.null))] = list(-1);
        conditional.probs.fit[unlist(lapply(conditional.probs.fit, is.null))] = list(1);

        ## Perform the estimation of the probabilities if requested.
        
        if (do.estimation) {
            ## Estimate the error rates and, given them, the
            ## probabilities for the causal topology.
            
            estimated.error.rates.fit =
                estimate.dag.error.rates(dataset,
                                         prima.facie.parents$marginal.probs,
                                         prima.facie.parents$joint.probs,
                                         parents.pos.fit);
            estimated.probabilities.fit =
                estimate.dag.probs(dataset,
                                   prima.facie.parents$marginal.probs,
                                   prima.facie.parents$joint.probs,
                                   parents.pos.fit,
                                   estimated.error.rates.fit);
        } else {
            estimated.error.rates.fit =
                list(error.fp = NA,
                     error.fn = NA);
            estimated.probabilities.fit =
                list(marginal.probs = NA,
                     joint.probs = NA,
                     conditional.probs = NA);
        }

        ## Set results for the current regolarizator.
        
        probabilities.observed =
            list(marginal.probs = prima.facie.parents$marginal.probs,
                 joint.probs = prima.facie.parents$joint.probs,
                 conditional.probs = conditional.probs.fit);
        probabilities.fit =
            list(estimated.marginal.probs = estimated.probabilities.fit$marginal.probs,
                 estimated.joint.probs = estimated.probabilities.fit$joint.probs,
                 estimated.conditional.probs = estimated.probabilities.fit$conditional.probs);
        probabilities =
            list(probabilities.observed = probabilities.observed,
                 probabilities.fit = probabilities.fit);
        parents.pos = parents.pos.fit;
        error.rates = estimated.error.rates.fit;

        ## Save the results for the model.
        
        model[[reg]] =
            list(probabilities = probabilities,
                 parents.pos = parents.pos,
                 error.rates = error.rates,
                 adj.matrix = best.parents$adj.matrix);
    }

    ## Set the execution parameters.
    
    parameters =
        list(algorithm = "CHOW_LIU",
             regularization = regularization,
             do.boot = do.boot,
             nboot = nboot,
             pvalue = pvalue,
             min.boot = min.boot,
             min.stat = min.stat,
             boot.seed = boot.seed,
             do.estimation = do.estimation,
             silent = silent);

    ## Return the results.
    
    topology =
        list(dataset = dataset,
             hypotheses = NA,
             adj.matrix.prima.facie = adj.matrix.prima.facie,
             confidence = prima.facie.parents$pf.confidence,
             model = model,
             parameters = parameters,
             execution.time = (proc.time() - ptm));
    return(topology);
}


# reconstruct the best causal topology by Chow Liu algorithm combined with probabilistic causation
# @title perform.likelihood.fit.chow.liu
# @param dataset a valid dataset
# @param adj.matrix the adjacency matrix of the prima facie causes
# @param regularization regularization term to be used in the likelihood fit
# @return topology: the adjacency matrix of both the prima facie and causal topologies
#
perform.likelihood.fit.chow.liu = function( dataset, adj.matrix, regularization ) {{

    ## Each variable should at least have 2 values: I'm ignoring
    ## connection to invalid events but, still, need to make the
    ## dataset valid for bnlearn.
    
    for (i in 1:ncol(dataset)) {
        if (sum(dataset[, i]) == 0) {
            dataset[sample(1:nrow(dataset), size=1), i] = 1;
        } else if (sum(dataset[, i]) == nrow(dataset)) {
            dataset[sample(1:nrow(dataset), size=1), i] = 0;
        }
    }

    # adjacency matrix of the topology reconstructed by likelihood fit
    adj.matrix.fit = array(0,c(nrow(adj.matrix),ncol(adj.matrix)))
    rownames(adj.matrix.fit) = colnames(dataset)
    colnames(adj.matrix.fit) = colnames(dataset)
    
    # set the regularizator
    if(regularization=="loglik") {
        regularization = "LR"
    }
    else if(regularization=="aic") {
        regularization = "AIC"
    }
    else if(regularization=="bic") {
        regularization = "BIC"
    }
    
    # set edges among disconnected nodes in the blacklist
    blacklist = NULL
    for (i in 1:nrow(adj.matrix)) {
        for (j in i:ncol(adj.matrix)) {
            if(i!=j && adj.matrix[i,j]==0 && adj.matrix[j,i]==0) {
                new_edge = c(i,j)
                blacklist = rbind(blacklist,new_edge)
            }
        }
    }
    
    # compute the best Chow-Liu tree among the valid edges
    best_chow_liu_tree = minForest(dataset,forbEdges=blacklist,stat=regularization)


    # get the best topology considering both the priors and the Chow-Liu tree
    if(length(best_chow_liu_tree@edges)>0) {
        for (i in 1:nrow(best_chow_liu_tree@edges)) {
            if(adj.matrix[best_chow_liu_tree@edges[i,1],best_chow_liu_tree@edges[i,2]]==1) {
                adj.matrix.fit[best_chow_liu_tree@edges[i,1],best_chow_liu_tree@edges[i,2]] = 1
            }
            else if(adj.matrix[best_chow_liu_tree@edges[i,2],best_chow_liu_tree@edges[i,1]]==1) {
                adj.matrix.fit[best_chow_liu_tree@edges[i,2],best_chow_liu_tree@edges[i,1]] = 1
            }
        }
    }
    
    ## Save the results and return them.
    
    adj.matrix =
        list(adj.matrix.pf = adj.matrix,
             adj.matrix.fit = adj.matrix.fit);
    topology = list(adj.matrix = adj.matrix);
    return(topology)

}


#### end of file -- chow.liu.algorithm.R
