#### TRONCO: a tool for TRanslational ONCOlogy
####
#### Copyright (c) 2015-2016, Marco Antoniotti, Giulio Caravagna, Luca De Sano,
#### Alex Graudenzi, Giancarlo Mauri, Bud Mishra and Daniele Ramazzotti.
####
#### All rights reserved. This program and the accompanying materials
#### are made available under the terms of the GNU GPL v3.0
#### which accompanies this distribution.

# reconstruct the best topology based on probabilistic causation and maximum likelihood estimation
# @title suppes.mle.fit
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
mltree.fit <- function(dataset,
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
            madonna.troia.wrapper(dataset,
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
        prima.facie.parents$adj.matrix$adj.matrix.acyclic

    matriciona = adj.matrix.prima.facie
    for (i in 1:nrow(matriciona)) {
        for (j in 1:ncol(matriciona)) {
            if (i != j) {
                matriciona[[i,j]] = 1
            }
        }
    }

    ## Perform the likelihood fit with the required strategy.
    
    model = list();
    for (reg in regularization) {

        ## Perform the likelihood fit with the chosen regularization
        ## score on the prima facie topology.
        
        if (!silent)
            cat('*** Performing likelihood-fit with regularization:', reg, '.\n')
        best.parents =
            compute.mltree(dataset,
                           adj.matrix.prima.facie,
                           regularization = reg)

        ## Set the structure to save the conditional probabilities of
        ## the reconstructed topology.

        reconstructed.model = create.model(dataset,
            best.parents,
            prima.facie.parents)

        model.name = paste('mltree', reg, sep='_')
        model[[model.name]] = reconstructed.model
    }

    ## Set the execution parameters.
    
    parameters =
        list(algorithm = "MLTREE",
             regularization = regularization,
             do.boot = do.boot,
             nboot = nboot,
             pvalue = pvalue,
             min.boot = min.boot,
             min.stat = min.stat,
             boot.seed = boot.seed,
             silent = silent);

    ## Return the results.
    
    topology =
        list(dataset = dataset,
             hypotheses = NA,
             adj.matrix.prima.facie = adj.matrix.prima.facie,
             confidence = prima.facie.parents$pf.confidence,
             model = model,
             tree.list = best.parents$tree.list,
             parameters = parameters,
             execution.time = (proc.time() - ptm))
    topology = rename.reconstruction.fields(topology, dataset)
    return(topology)
}


# reconstruct the best causal topology by maximum likelihood estimation combined with probabilistic causation
# @title perform.likelihood.fit.mle
# @param dataset a valid dataset
# @param adj.matrix the adjacency matrix of the prima facie causes
# @param regularization regularization term to be used in the likelihood fit
# @param command type of search, either hill climbing (hc) or tabu (tabu)
# @return topology: the adjacency matrix of both the prima facie and causal topologies
#
compute.mltree = function(dataset,
                          adj.matrix,
                          regularization) {

    data = as.categorical.dataset(dataset)
    adj.matrix.prima.facie = adj.matrix
    
    # adjacency matrix of the topology reconstructed by likelihood fit
    adj.matrix.fit = array(0,c(nrow(adj.matrix),ncol(adj.matrix)))
    rownames(adj.matrix.fit) = colnames(dataset)
    colnames(adj.matrix.fit) = colnames(dataset)
    
    internal.which <- function(x) {
        if (sum(adj.matrix.prima.facie[,x]) > 0) 
            return(paste0('which(adj.matrix.prima.facie[,', x,'] == 1)'))
        else 
            return('0')
    }

    com = paste(sapply(1:nrow(adj.matrix.prima.facie), internal.which))
    com = paste0(com, collapse= ', ')
    com = paste0('expand.grid(', com,')', collapse = '')

    indexes = eval(parse(text=com))


    #print(indexes)

    # Generate all matrices - amazing
    candidates = Reduce(
        append,
        apply(indexes, 1, function(x) {
            M = matrix(0, nrow=ncol(adj.matrix.prima.facie), ncol = ncol(adj.matrix.prima.facie))
            rownames(M) = colnames(dataset)
            colnames(M) = colnames(dataset)
            for(i in 1:nrow(M)) {
                M[x[i], i] = 1
            }
            #print(M)
            M.graph = graph.adjacency(M)
            #print(is.dag(M.graph))
            if (is.dag(M.graph)) {
                return(list(M))
            } else {
                return(list())
            }
        })
    )
    cat("Total number of prima facie trees: ", length(candidates), '\n')

    MAXRESULTS = 50 # output results -- full posterior by setting it to NUMTESTS
    NUMTESTS = length(candidates)

    tree.list = NULL

    if (NUMTESTS > 1) {
        # empty networks
        net = empty.graph(colnames(dataset), num = NUMTESTS)

        # set edges/compute score
        scores = rep(0, NUMTESTS)
        # make this parallel

        for(i in 1:NUMTESTS) {
            amat(net[[i]]) = candidates[[i]]
            scores[i] = score(net[[i]], data, type = "loglik")
        }

        # sort and trim results
        net = net[order(scores, decreasing = TRUE)]
        net = net[1: min(MAXRESULTS, length(net))]

        best.net = net[[1]]
        adj.matrix.fit = amat(best.net)
        for (i in 1:length(net)) {
            tree.list[[i]] = amat(net[[i]])
        }


    } else {
        adj.matrix.fit = candidates[[1]]
        tree.list[[1]] = candidates[[1]]
    }
    
    ## Save the results and return them.
    adj.matrix =
        list(adj.matrix.pf = adj.matrix.prima.facie,
             adj.matrix.fit = adj.matrix.fit)
    topology = list(adj.matrix = adj.matrix, tree.list = tree.list)

    return(topology)
}

#### end of file -- mltree.algorithm.R



madonna.troia.wrapper <- function(dataset,
                                            hypotheses,
                                            nboot,
                                            pvalue,
                                            adj.matrix,
                                            min.boot,
                                            min.stat,
                                            boot.seed,
                                            silent ) {

    ## Perform a robust estimation of the scores using rejection
    ## sampling bootstrap.
    
    scores =
        get.bootstrapped.scores(dataset,
                               nboot,adj.matrix,
                               min.boot,
                               min.stat,
                               boot.seed,
                               silent,
                               0.0,
                               0.0);

    ## Compute the observed and joint probabilities as the mean of the
    ## bootstrapped values.
    
    marginal.probs = array(-1, dim = c(ncol(dataset), 1));
    joint.probs = array(-1, dim = c(ncol(dataset), ncol(dataset)));
    for (i in 1:ncol(dataset)) {
        marginal.probs[i, 1] =
            mean(unlist(scores$marginal.probs.distributions[i, 1]));
        for (j in i:ncol(dataset)) {
            joint.probs[i,j] =
                mean(unlist(scores$joint.probs.distributions[i, j]));
            if (i != j) {
                joint.probs[j, i] = joint.probs[i, j];
            }
        }
    }

    ## Remove all the edges not representing a prima facie cause.
    
    prima.facie.topology =
        madonna.troia(adj.matrix,
                                       hypotheses,
                                       scores$marginal.probs.distributions,
                                       scores$prima.facie.model.distributions,
                                       scores$prima.facie.null.distributions,
                                       pvalue,
                                       dataset,
                                       marginal.probs,joint.probs,
                                       silent);

    ## Save the results and return them.
    
    prima.facie.parents <-
        list(marginal.probs = marginal.probs,
             joint.probs = joint.probs,
             adj.matrix = prima.facie.topology$adj.matrix,
             pf.confidence = prima.facie.topology$edge.confidence.matrix);
    return(prima.facie.parents);
}






# select the best set of prima facie causes per node
# @title get.prima.facie.causes.do.boot
# @param adj.matrix adjacency matrix of the initially valid edges
# @param hypotheses hypotheses to be considered
# @param marginal.probs.distributions distributions of the bootstrapped marginal probabilities
# @param prima.facie.model.distributions distributions of the prima facie model
# @param prima.facie.null.distributions distributions of the prima facie null
# @param pvalue minimum pvalue for the Mann-Whitney U tests to be significant
# @param dataset a valid dataset
# @param marginal.probs observed marginal probabilities
# @param joint.probs observed joint probabilities
# @param silent Should I be verbose?
# @return prima.facie.topology: list describing the topology of the prima facie causes
#
madonna.troia <- function(adj.matrix,
                                           hypotheses,
                                           marginal.probs.distributions,
                                           prima.facie.model.distributions,
                                           prima.facie.null.distributions,
                                           pvalue,
                                           dataset,
                                           marginal.probs,
                                           joint.probs,
                                           silent = FALSE) {

    ## Structure to save the confidence of the edges.
    
    edge.confidence.matrix <- array(list(), c(3, 1));
    edge.confidence.matrix[[1,1]] =
        array(NA, c(ncol(prima.facie.model.distributions),
                    ncol(prima.facie.model.distributions)));
    edge.confidence.matrix[[2,1]] =
        array(NA, c(ncol(prima.facie.model.distributions),
                    ncol(prima.facie.model.distributions)));
    edge.confidence.matrix[[3,1]] =
        array(NA, c(ncol(prima.facie.model.distributions),
                    ncol(prima.facie.model.distributions)));

    ## Verify Suppes' conditions for prima facie causes;
    ## i.e., i --> j implies P(i)>P(j) (temporal priority) and
    ## P(j|i)>P(j|not i) (probability raising).
    
    ## Verify the temporal priority condition.
    
    if (!silent)
        cat(paste0('\tEvaluating \"temporal priority\" (Wilcoxon, p-value ',
                   pvalue,
                   ')\n'));
    temporal.priority =
        verify.temporal.priority.do.boot(marginal.probs.distributions,
                                         pvalue,adj.matrix,
                                         edge.confidence.matrix);

    ## Verify the probability raising condition.
    
    if (!silent)
        cat(paste0('\tEvaluating \"probability raising\" (Wilcoxon, p-value ',
                   pvalue,
                   ')\n'));
    probability.raising =
        verify.probability.raising.do.boot(prima.facie.model.distributions,
                                           prima.facie.null.distributions,
                                           pvalue,temporal.priority$adj.matrix,temporal.priority$edge.confidence.matrix);

    ## Perform the hypergeometric test for each pair of events.
    
    for (i in 1:ncol(adj.matrix)) {
        for (j in i:nrow(adj.matrix)) {

            ## The diagonal (self cause) and the other invalid edges
            ## have not to be considered.
            
            if (adj.matrix[i, j] != 0 || adj.matrix[j, i] != 0) {
                
                ## Compute the confidence by hypergeometric test for
                ## both j --> i and i --> j.
                
                probability.raising$edge.confidence.matrix[[3, 1]][i, j] =
                    phyper(joint.probs[i, j] * nrow(dataset),
                           marginal.probs[i] * nrow(dataset),
                           nrow(dataset) - marginal.probs[i] * nrow(dataset),
                           marginal.probs[j] * nrow(dataset),
                           lower.tail = FALSE);
                probability.raising$edge.confidence.matrix[[3, 1]][j, i] =
                    probability.raising$edge.confidence.matrix[[3, 1]][i, j];
            } else {
                probability.raising$edge.confidence.matrix[[3, 1]][i, j] = 1;
                probability.raising$edge.confidence.matrix[[3, 1]][j, i] = 1;
            }
        }
    }

    ## Remove any cycle.
    
    #adj.matrix.cyclic = probability.raising$adj.matrix
    adj.matrix.cyclic = temporal.priority$adj.matrix
#    if (length(temporal.priority$not.ordered) > 0
#        || !is.na(hypotheses[1])) {
#
#        if (!silent)
#            cat('*** Loop detection found loops to break.\n')
#
#        weights.temporal.priority =
#            probability.raising$edge.confidence.matrix[[1, 1]] +
#                probability.raising$edge.confidence.matrix[[2, 1]];
#        weights.matrix =
#            probability.raising$edge.confidence.matrix[[2, 1]] +
#                probability.raising$edge.confidence.matrix[[3, 1]];
#        acyclic.topology =
#            remove.cycles(probability.raising$adj.matrix,
#                          weights.temporal.priority,
#                          weights.matrix,
#                          temporal.priority$not.ordered,
#                          hypotheses,
#                          silent);
#        adj.matrix.acyclic = acyclic.topology$adj.matrix;
#
#    } else {
#        adj.matrix.acyclic = probability.raising$adj.matrix;
#    }

    adj.matrix.acyclic = temporal.priority$adj.matrix

    adj.matrix =
        list(adj.matrix.cyclic = adj.matrix.cyclic,
             adj.matrix.acyclic = adj.matrix.acyclic)

    ## Save the results and return them.
    
    prima.facie.topology =
        list(adj.matrix = adj.matrix,
             edge.confidence.matrix = probability.raising$edge.confidence.matrix);
    return(prima.facie.topology);
}
