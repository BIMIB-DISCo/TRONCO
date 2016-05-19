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
        prima.facie.parents$adj.matrix$adj.matrix.acyclic

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


    # Generate all matrices - amazing
    candidates = Reduce(
        append,
        apply(indexes, 1, function(x) {
            M = matrix(0L, nrow=ncol(adj.matrix.prima.facie), ncol = ncol(adj.matrix.prima.facie))
            rownames(M) = colnames(dataset)
            colnames(M) = colnames(dataset)
            for(i in 1:nrow(M)) {
                M[x[i], i] = 1L
            }
            return(list(M))
        })
    )
    cat("Total number of prima facie trees: ", length(candidates), '\n')

    MAXRESULTS = 10 # output results -- full posterior by setting it to NUMTESTS
    NUMTESTS = length(candidates)

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
        net = net[order(scores)]
        net = net[1: min(MAXRESULTS, length(net))]

        best.net = net[[1]]
        adj.matrix.fit = amat(best.net)
    } else {
        adj.matrix.fit = candidates[[1]]
    }
    
    ## Save the results and return them.
    adj.matrix =
        list(adj.matrix.pf = adj.matrix.prima.facie,
             adj.matrix.fit = adj.matrix.fit)
    topology = list(adj.matrix = adj.matrix)

    return(topology)
}

#### end of file -- mle.algorithm.R
