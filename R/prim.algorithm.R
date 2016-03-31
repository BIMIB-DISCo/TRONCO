#### TRONCO: a tool for TRanslational ONCOlogy
####
#### Copyright (c) 2015-2016, Marco Antoniotti, Giulio Caravagna, Luca De Sano,
#### Alex Graudenzi, Giancarlo Mauri, Bud Mishra and Daniele Ramazzotti.
####
#### All rights reserved. This program and the accompanying materials
#### are made available under the terms of the GNU GPL v3.0
#### which accompanies this distribution.


# reconstruct the best topology based on probabilistic causation and Prim algorithm
# @title prim.fit
# @param dataset a dataset describing a progressive phenomenon
# @param regularization regularizators to be used for the likelihood fit
# @param do.boot should I perform bootstrap? Yes if TRUE, no otherwise
# @param nboot integer number (greater than 0) of bootstrap sampling to be performed
# @param pvalue pvalue for the tests (value between 0 and 1)
# @param min.boot minimum number of bootstrapping to be performed
# @param min.stat should I keep bootstrapping untill I have nboot valid values?
# @param boot.seed seed to be used for the sampling
# @param silent should I be verbose?
# @return topology: the reconstructed tree topology
#
prim.fit <- function(dataset,
                      regularization = "no_reg",
                      do.boot = TRUE,
                      nboot = 100,
                      pvalue = 0.05,
                      min.boot = 3,
                      min.stat = TRUE,
                      boot.seed = NULL,
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
            cat('*** Performing likelihood-fit with regularization:',reg,'.\n')
        best.parents =
            perform.likelihood.fit.prim(dataset,
                                   prima.facie.parents$adj.matrix$adj.matrix.acyclic,
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
                        list(c(unlist(parents.pos.fit[j, 1]), i))
                    conditional.probs.fit[j,1] =
                        list(c(unlist(conditional.probs.fit[j, 1]),
                               prima.facie.parents$joint.probs[i, j] /
                                   prima.facie.parents$marginal.probs[i]))
                }
            }
        }
        parents.pos.fit[unlist(lapply(parents.pos.fit, is.null))] = list(-1)
        conditional.probs.fit[unlist(lapply(conditional.probs.fit, is.null))] = list(1)

        ## Perform the estimation of the probabilities if requested.
        
        estimated.error.rates.fit =
            list(error.fp = NA,
                 error.fn = NA)
        estimated.probabilities.fit =
            list(marginal.probs = NA,
                 joint.probs = NA,
                 conditional.probs = NA)

        ## Set results for the current regolarizator.
        
        probabilities.observed =
            list(marginal.probs = prima.facie.parents$marginal.probs,
                 joint.probs = prima.facie.parents$joint.probs,
                 conditional.probs = conditional.probs.fit)
        probabilities.fit =
            list(estimated.marginal.probs = estimated.probabilities.fit$marginal.probs,
                 estimated.joint.probs = estimated.probabilities.fit$joint.probs,
                 estimated.conditional.probs = estimated.probabilities.fit$conditional.probs)
        probabilities =
            list(probabilities.observed = probabilities.observed,
                 probabilities.fit = probabilities.fit)
        parents.pos = parents.pos.fit
        error.rates = estimated.error.rates.fit

        ## Save the results for the model.

        model.name = paste('prim', reg, sep='_')
        
        model[[model.name]] =
            list(probabilities = probabilities,
                 parents.pos = parents.pos,
                 error.rates = error.rates,
                 adj.matrix = best.parents$adj.matrix)
    }

    ## Set the execution parameters.
    
    parameters =
        list(algorithm = "PRIM",
             regularization = regularization,
             do.boot = do.boot,
             nboot = nboot,
             pvalue = pvalue,
             min.boot = min.boot,
             min.stat = min.stat,
             boot.seed = boot.seed,
             silent = silent)

    ## Return the results.
    
    topology =
        list(dataset = dataset,
             hypotheses = NA,
             adj.matrix.prima.facie = adj.matrix.prima.facie,
             confidence = prima.facie.parents$pf.confidence,
             model = model,
             parameters = parameters,
             execution.time = (proc.time() - ptm))
    topology = rename.reconstruction.fields(topology, dataset)
    return(topology)
}


# reconstruct the best causal topology by Prim algorithm combined with probabilistic causation
# title perform.likelihood.fit.prim
# param dataset a valid dataset
# param adj.matrix the adjacency matrix of the prima facie causes
# param regularization regularization term to be used in the likelihood fit
# param command type of search, either hill climbing (hc) or tabu (tabu)
# return topology: the adjacency matrix of both the prima facie and causal topologies
#
perform.likelihood.fit.prim = function(dataset,
                                       adj.matrix,
                                       regularization,
                                       command = "hc"){

    data = as.categorical.dataset(dataset)
    adj.matrix.prima.facie = adj.matrix
    
    # adjacency matrix of the topology reconstructed by likelihood fit
    adj.matrix.fit = array(0,c(nrow(adj.matrix),ncol(adj.matrix)))
    rownames(adj.matrix.fit) = colnames(dataset)
    colnames(adj.matrix.fit) = colnames(dataset)
       
    # create the graph of valid edges
    curr.valid.adj.matrix = array(0, c(nrow(adj.matrix), ncol(adj.matrix)))
    cont = 0
    for (i in 1:nrow(curr.valid.adj.matrix)) {
        for (j in i:nrow(curr.valid.adj.matrix)) {
            if (!(adj.matrix[i,j] == 0 && adj.matrix[j,i] == 0)) {
                cont = cont + 1
                curr.valid.adj.matrix[i,j] = 1
                curr.valid.adj.matrix[j,i] = 1
            }
        }
    }
    
    # if the graph has at least one edge, build to weighted igraph object
    if (cont > 0) {
        curr.graph = graph.adjacency(curr.valid.adj.matrix, mode = "undirected")
        all_edges = get.edgelist(curr.graph)
        # set the weights to the edges
        new_weights = NULL
        if (cont == 1) {
            new_weights = mutinformation(data[ ,all_edges[1,1]], data[ ,all_edges[1,2]])
        } else {
            for (i in 1:nrow(all_edges)) {
                new_weights = c(new_weights, 
                    mutinformation(data[ ,all_edges[i,1]], data[ ,all_edges[i,2]]))
            }
        }
        # set the weights to the graph
        E(curr.graph)$weight = max(new_weights) - new_weights
        # get the minimum spanning tree by Prim algorithm
        curr.valid.adj.matrix = as.matrix(get.adjacency(minimum.spanning.tree(curr.graph, 
                                                                              algorithm="prim")))

        # build the matrix of the priors
        new_prior_matrix = array(0, c(nrow(adj.matrix), ncol(adj.matrix)))
        rownames(new_prior_matrix) = colnames(dataset)
        colnames(new_prior_matrix) = colnames(dataset)
        for (i in 1:nrow(curr.valid.adj.matrix)) {
            for (j in i:ncol(curr.valid.adj.matrix)) {
                if (i != j) {
                    if (curr.valid.adj.matrix[i,j] == 1) {
                        if (adj.matrix[i,j] == 1) {
                            new_prior_matrix[i,j] = 1
                        } else if(adj.matrix[j,i] == 1) {
                            new_prior_matrix[j,i] = 1
                        }
                    }
                }
            }
        }
        adj.matrix = new_prior_matrix

        ## Perform the likelihood fit if requested
        adj.matrix.fit = lregfit(data,
            adj.matrix,
            adj.matrix.fit,
            regularization,
            command)
    }
    
    
    adj.matrix = list(adj.matrix.pf = adj.matrix.prima.facie,
        adj.matrix.fit = adj.matrix.fit)
    topology = list(adj.matrix = adj.matrix)
    return(topology)
}


#### end of file -- prim.algorithm.R
