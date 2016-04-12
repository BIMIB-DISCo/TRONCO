#### TRONCO: a tool for TRanslational ONCOlogy
####
#### Copyright (c) 2015-2016, Marco Antoniotti, Giulio Caravagna, Luca De Sano,
#### Alex Graudenzi, Giancarlo Mauri, Bud Mishra and Daniele Ramazzotti.
####
#### All rights reserved. This program and the accompanying materials
#### are made available under the terms of the GNU GPL v3.0
#### which accompanies this distribution.

# reconstruct the best dag topology running CAPRI algorithm
# @title capri.fit
# @param dataset a dataset describing a progressive phenomenon
# @param hypotheses hypotheses to be considered in the reconstruction
# @param command type of search for the likelihood fit, either hill climbing (hc) or tabu (tabu)
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
capri.fit <- function(dataset,
                      hypotheses = NA,
                      command = "hc",
                      regularization = c("bic", "aic"),
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

    ## Consider any hypothesis.
    
    adj.matrix = hypothesis.adj.matrix(hypotheses, adj.matrix);

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
                                            hypotheses,
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
                                            hypotheses,
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

    ## Perform the likelihood fit with the required regularization
    ## scores.
    
    model = list();
    for (reg in regularization) {

        ## Perform the likelihood fit with the chosen regularization
        ## score on the prima facie topology.
        
        if (!silent)
            cat(paste0('*** Performing likelihood-fit with regularization ',reg,'.\n'))
        best.parents =
            perform.likelihood.fit.capri(dataset,
                prima.facie.parents$adj.matrix$adj.matrix.acyclic,
                command,
                regularization = reg)

        ## Set the structure to save the conditional probabilities of
        ## the reconstructed topology.

        reconstructed.model = create.model(dataset,
            best.parents,
            prima.facie.parents)

        model.name = paste('capri', reg, sep='_')
        model[[model.name]] = reconstructed.model
    }

    ## Set the execution parameters.
    
    parameters =
        list(algorithm = "CAPRI",
             command = command,
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
             hypotheses = hypotheses,
             adj.matrix.prima.facie = adj.matrix.prima.facie,
             confidence = prima.facie.parents$pf.confidence,
             model = model,
             parameters = parameters,
             execution.time = (proc.time() - ptm))
    topology = rename.reconstruction.fields(topology, dataset)
    return(topology)
}


# check if the dataset is valid accordingly to the probability raising
# @title check.dataset
# @param dataset a dataset describing a progressive phenomenon
# @param adj.matrix adjacency matrix of the topology
# @param verbose should I print the warnings? Yes if TRUE, no otherwise
# @return valid.dataset: a dataset valid accordingly to the probability raising
check.dataset <- function(dataset, adj.matrix, verbose ) {

    ## Perform the preprocessing only if I have at least two binary
    ## events and two samples.
    
    if (length(ncol(dataset)) > 0
        && ncol(dataset) > 1
        && length(nrow(dataset)) > 0
        && nrow(dataset) > 1
        && length(dataset[dataset == 0 | dataset == 1]) == nrow(dataset) * ncol(dataset)) {

        ## Structure to compute the observed and observed joint
        ## probabilities.
        
        pair.count <- array(0, dim = c(ncol(dataset), ncol(dataset)));
        
        ## Compute the probabilities on the dataset.
        
        for (i in 1:ncol(dataset)) {
            for (j in 1:ncol(dataset)) {
                val1 = dataset[ ,i];
                val2 = dataset[ ,j];
                pair.count[i,j] = (t(val1) %*% val2);
            }
        }
        ## Marginal.probs is an array of the observed marginal
        ## probabilities.
        
        marginal.probs =
            array(as.matrix(diag(pair.count) / nrow(dataset)),
                  dim = c(ncol(dataset),1));
        
        ## joint.probs is an array of the observed joint
        ## probabilities.
        
        joint.probs = as.matrix(pair.count/nrow(dataset));

        ## Evaluate the connections.
        
        invalid.events = vector();
        for (i in 1:ncol(adj.matrix)) {
            for (j in 1:nrow(adj.matrix)) {

                ## if i --> j is valid
                if (i != j && adj.matrix[i, j] == 1) {
                    
                    ## the potential cause is always present
                    
                    if (marginal.probs[i] == 1) {
                        
                        ## the potential child is not always missing
                        
                        if (marginal.probs[i] > 0) {
                            adj.matrix[i, j] = 0;
                            ## invalid.events = rbind(invalid.events,t(c(i,j)));
                        }
                        ## the potential child is always missing
                        else if (marginal.probs[i] == 0) {
                            adj.matrix[i, j] = 0;
                        }
                    }
                    ## the potential cause is always missing
                    else if (marginal.probs[i] == 0) {
                        adj.matrix[i, j] = 0;
                    }
                    ## the potential child is always present
                    else if (marginal.probs[j] == 1) {
                        adj.matrix[i, j] = 0;
                    }
                    ## the potential child is always missing
                    else if (marginal.probs[j] == 0) {
                        adj.matrix[i, j] = 0;
                    }
                    ## the two events are equals
                    else if ((joint.probs[i,j] / marginal.probs[i]) == 1
                             && (joint.probs[i,j] / marginal.probs[j]) == 1) {
                        adj.matrix[i, j] = 0;
                        invalid.events = rbind(invalid.events,t(c(i,j)));
                    }
                }
            }
        }
        
        if (length(invalid.events) > 0) {
            colnames(invalid.events) = c("cause","effect");
        }
        valid.dataset =
            list(dataset = dataset,
                 adj.matrix = adj.matrix,
                 invalid.events = invalid.events,
                 marginal.probs = marginal.probs,
                 joint.probs = joint.probs);
    }

    ## if the dataset is not valid, we stop here
    else {
        if (verbose == TRUE) {
            warning("The dataset must contain at least two binary events and two samples.");
        }
        valid.dataset =
            list(dataset = NA,
                 adj.matrix = NA,
                 invalid.events = NA,
                 marginal.probs = NA,
                 joint.probs = NA);
    }
    return(valid.dataset);
}


# compute a robust estimation of the scores using rejection sampling bootstrap
# @title get.bootstapped.scores
# @param dataset a valid dataset
# @param nboot number of bootstrap resampling to be performed
# @param adj.matrix adjacency matrix of the initially valid edges
# @param min.boot minimum number of bootstrapping to be performed
# @param min.stat should I keep bootstrapping untill I have nboot valid values?
# @param boot.seed seed to be used for the sampling
# @param silent Should I be verbose?
# @return scores: list structure with the scores and the data generated by bootstrap
#
get.bootstapped.scores <- function(dataset,
                                   nboot,
                                   adj.matrix,
                                   min.boot = 3,
                                   min.stat = TRUE,
                                   boot.seed = NULL,
                                   silent = FALSE) {

    ## Structures to save the distributions generated by the
    ## bootstrapped datasets.
    
    marginal.probs.distributions <-
        array(list(-1), c(ncol(dataset), 1));
    joint.probs.distributions <-
        array(list(-1), c(ncol(dataset), ncol(dataset)));
    prima.facie.model.distributions <-
        array(list(-1), c(ncol(dataset), ncol(dataset)));
    prima.facie.null.distributions <-
        array(list(-1), c(ncol(dataset), ncol(dataset)));
    
    ## Structures to save the number of performed valid (not rejected)
    ## sampling.
    
    sampled.marginal.probs.distributions <-
        array(0, dim = c(ncol(dataset),1));
    sampled.joint.probs.distributions <-
        array(0, dim = c(ncol(dataset), ncol(dataset)));
    sampled.prima.facie.distributions <-
        array(0, dim = c(ncol(dataset),ncol(dataset)));

    ## I require a minimum of min.boot (default = 3) sampling of
    ## bootstrap.
    
    nboot = max(nboot,min.boot);

    ## Set not to sample for the invalid edges.
    
    for (i in 1:nrow(adj.matrix)) {
        for (j in 1:ncol(adj.matrix)) {
            if (adj.matrix[i, j] == 0) {
                sampled.prima.facie.distributions[i,j] = nboot;
            }
        }
    }

    ## Perform bootstrap estimation based on a number of bootstrapped
    ## (>= nboot) datasets.
    
    curr.iteration = min(sampled.prima.facie.distributions);
    boot.counter = 0;

    ## Set the seed to be used for the sampling.
    
    set.seed(boot.seed);

    #if (silent == FALSE) {
    #    ## Create a progress bar.
    #    flush.console();
    #    pb <- txtProgressBar(curr.iteration, nboot, style = 3);
    #}

    dot = 0
    if (!silent) {
        cat('\t')
    }

    while (curr.iteration<nboot) {

        ## Define the dataset to be used in this iteration and compute
        ## the scores on it.
        
        sampled.data =
            dataset[sample(1:nrow(dataset),
                           size = nrow(dataset),
                           replace = TRUE),
                    ];
        boot.counter = boot.counter + 1;

        ## Compute the scores on the sampled data.
        
        curr.scores = get.dag.scores(sampled.data,adj.matrix);
        curr.marginal.probs = curr.scores$marginal.probs;
        curr.joint.probs = curr.scores$joint.probs;
        curr.prima.facie.model = curr.scores$prima.facie.model;
        curr.prima.facie.null = curr.scores$prima.facie.null;

        ## Save the (valid) scores for each edge.
        
        for (i in 1:nrow(curr.prima.facie.model)) {

            ## Get the marginal probabilities from the sampled data.
            
            if (sampled.marginal.probs.distributions[i, 1] == 0) {
                sampled.marginal.probs.distributions[i, 1] = 1;
                marginal.probs.distributions[i,1] =
                    curr.marginal.probs[i, 1];
            } else {
                marginal.probs.distributions[i, 1] =
                    list(c(unlist(marginal.probs.distributions[i, 1]),
                           curr.marginal.probs[i, 1]));
            }

            for (j in 1:ncol(curr.prima.facie.model)) {

                ## Get the joint probs from the sampled data.
                
                if (sampled.joint.probs.distributions[i, j] == 0) {
                    sampled.joint.probs.distributions[i, j] = 1;
                    joint.probs.distributions[i, j] =
                        curr.joint.probs[i, j];
                } else {
                    joint.probs.distributions[i, j] =
                        list(c(unlist(joint.probs.distributions[i, j]),
                               curr.joint.probs[i, j]));
                }

                ## Get the prima facie estimations from the sampled
                ## data.
                
                if (curr.prima.facie.model[i, j] != -1) {
                    
                    ## Count the valid values per edge.
                    
                    sampled.prima.facie.distributions[i, j] =
                        sampled.prima.facie.distributions[i, j] + 1;

                    ## Save the scores.
                    
                    if (sampled.prima.facie.distributions[i, j] == 1) {
                        ## scores for i --> j
                        prima.facie.model.distributions[i, j] =
                            list(curr.prima.facie.model[i, j]);
                        prima.facie.null.distributions[i, j] =
                            list(curr.prima.facie.null[i, j]);
                    } else {
                        ## Scores for i --> j
                        prima.facie.model.distributions[i, j] =
                            list(c(unlist(prima.facie.model.distributions[i, j]),
                                   curr.prima.facie.model[i, j]));
                        prima.facie.null.distributions[i, j] =
                            list(c(unlist(prima.facie.null.distributions[i, j]),
                                   curr.prima.facie.null[i, j]));
                    }
                }
            }
        }

        ## Set the number of performed iterations after the last
        ## bootstrap sampling.
        
        curr.iteration = min(sampled.prima.facie.distributions);

        ## If the flag min.stat is FALSE, even if after nboot
        ## iterations I don't have nboot valid entries, as soon as I
        ## have at least min.boot values, I stop anyway.
        
        if (min.stat == FALSE
            && boot.counter >= nboot
            && curr.iteration >= min.boot) {
            curr.iteration = nboot;
        }

        #if (silent == FALSE) {
        #    ## Increment the progress bar.
        #    if (min.stat == FALSE) {
        #        setTxtProgressBar(pb, boot.counter);
        #    } else {
        #        setTxtProgressBar(pb, curr.iteration);
        #    }
        #}
        if (!silent && (curr.iteration %% 5 == 0)) {
            cat('.')
            dot = dot + 1
            if (dot %% 50 == 0) {
                cat('\n\t')
            }
        }
    }

    #if (silent == FALSE) {
    #    ## Close the progress bar.
    #    close(pb);
    #}
    if (!silent) {
        cat('\n')
    }

    ## Save the results and return them.
    scores =
        list(marginal.probs.distributions = marginal.probs.distributions,
             joint.probs.distributions = joint.probs.distributions,
             prima.facie.model.distributions = prima.facie.model.distributions,
             prima.facie.null.distributions = prima.facie.null.distributions);
    return(scores);
}


# compute the observed probabilities and the prima facie scores on the dataset
# @title get.dag.scores
# @param dataset a valid dataset
# @param adj.matrix adjacency matrix of the initially valid edges
# @return scores: observed probabilities and prima facie scores
#
get.dag.scores <- function( dataset, adj.matrix ) {

    ## Structure to save the prima facie scores.
    
    prima.facie.model <- array(-1, dim = c(ncol(dataset), ncol(dataset)));
    prima.facie.null  <- array(-1, dim = c(ncol(dataset), ncol(dataset)));

    ## Structure to save the observed and observed-joint
    ## probabilities.
    
    pair.count <- array(0, dim = c(ncol(dataset), ncol(dataset)));

    ## Compute the observed probabilities on the dataset.
    
    for (i in 1:ncol(dataset)) {
        for (j in 1:ncol(dataset)) {
            val1 = dataset[ ,i];
            val2 = dataset[ ,j];
            pair.count[i,j] = (t(val1) %*% val2);
        }
    }
    
    ## marginal.probs is an array with the marginal probabilities.
    
    marginal.probs <-
        array(as.matrix(diag(pair.count) / nrow(dataset)),
              dim = c(ncol(dataset), 1));
    
    ## joint.probs is an array with the joint observed probabilities.
    
    joint.probs <- as.matrix(pair.count / nrow(dataset));

    ## Compute the prima facie scores based on the probability raising
    ## model.
    
    for (i in 1:nrow(prima.facie.model)) {
        for (j in 1:ncol(prima.facie.model)) {
            
            ## The scores are saved in the convention of the adjacency
            ## matrix, i.e., [i,j] means i is causing j  the diagonal
            ## (self cause) and the other invalid edges have not to be
            ## considered.
            
            if (adj.matrix[i, j] != 0) {
                
                ## Check if the connections from j to i and from i to
                ## j can be evaluated on this dataset.
                
                if (marginal.probs[i] > 0
                    && marginal.probs[i] < 1
                    && marginal.probs[j] > 0
                    && marginal.probs[j] < 1) {
                    
                    ## Check if the two events i and j are
                    ## distinguishable.
                    
                    if ((joint.probs[i, j] / marginal.probs[j]) < 1
                        || (joint.probs[i, j] / marginal.probs[i]) < 1) {
                        
                        ## prima facie scores of i --> j
                        
                        prima.facie.model[i, j] =
                            joint.probs[j, i] / marginal.probs[i];
                        prima.facie.null[i, j] =
                            (marginal.probs[j] - joint.probs[j, i]) / (1 - marginal.probs[i]);
                    }
                }
            }
        }
    }

    ## Save the results and return them.
    scores =
        list(marginal.probs = marginal.probs,
             joint.probs = joint.probs,
             prima.facie.model = prima.facie.model,
             prima.facie.null = prima.facie.null);
    return(scores);
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
get.prima.facie.causes.do.boot <- function(adj.matrix,
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
    
    adj.matrix.cyclic = probability.raising$adj.matrix
    if (length(temporal.priority$not.ordered) > 0
        || !is.na(hypotheses[1])) {

        if (!silent)
            cat('*** Loop detection found loops to break.\n')

        weights.temporal.priority =
            probability.raising$edge.confidence.matrix[[1, 1]] +
                probability.raising$edge.confidence.matrix[[2, 1]];
        weights.matrix =
            probability.raising$edge.confidence.matrix[[2, 1]] +
                probability.raising$edge.confidence.matrix[[3, 1]];
        acyclic.topology =
            remove.cycles(probability.raising$adj.matrix,
                          weights.temporal.priority,
                          weights.matrix,
                          temporal.priority$not.ordered,
                          hypotheses,
                          silent);
        adj.matrix.acyclic = acyclic.topology$adj.matrix;

    } else {
        adj.matrix.acyclic = probability.raising$adj.matrix;
    }
    adj.matrix =
        list(adj.matrix.cyclic = adj.matrix.cyclic,
             adj.matrix.acyclic = adj.matrix.acyclic)

    ## Save the results and return them.
    
    prima.facie.topology =
        list(adj.matrix = adj.matrix,
             edge.confidence.matrix = probability.raising$edge.confidence.matrix);
    return(prima.facie.topology);
}


# select the best set of prima facie causes per node without bootstrap
# @title get.prima.facie.causes.no.boot
# @param adj.matrix adjacency matrix of the initially valid edges
# @param hypotheses hypotheses object related to adjacency matrix
# @param marginal.probs observed marginal probabilities
# @param prima.facie.model prima facie model
# @param prima.facie.null prima facie null
# @param dataset a valid dataset
# @param joint.probs observed joint probabilities
# @param silent Should I be verbose?
# @return prima.facie.topology: adjacency matrix of the prima facie causes
#
get.prima.facie.causes.no.boot <- function(adj.matrix,
                                           hypotheses,
                                           marginal.probs,
                                           prima.facie.model,
                                           prima.facie.null,
                                           dataset,
                                           joint.probs,
                                           silent = FALSE) {

    ## Structure to save the confidence of the edges.
    
    edge.confidence.matrix <- array(list(), c(3, 1));
    edge.confidence.matrix[[1, 1]] = array(NA, c(ncol(adj.matrix), ncol(adj.matrix)));
    edge.confidence.matrix[[2, 1]] = array(NA, c(ncol(adj.matrix), ncol(adj.matrix)));
    edge.confidence.matrix[[3, 1]] = array(NA, c(ncol(adj.matrix), ncol(adj.matrix)));

    ## Verify Suppes' conditions for prima facie causes;
    ## i.e., i --> j implies P(i)>P(j) (temporal priority) and
    ## P(j|i)>P(j|not i) (probability raising).
    
    ## Verify the temporal priority condition.
    
    if (!silent)
        cat(paste0('\tEvaluating \"temporal priority\".\n'));
    temporal.priority =
        verify.temporal.priority.no.boot(marginal.probs,
                                         adj.matrix,
                                         edge.confidence.matrix);

    ## Verify the probability raising condition.
    
    if (!silent)
        cat(paste0('\tEvaluating \"probability raising\".\n'));
    probability.raising =
        verify.probability.raising.no.boot(prima.facie.model,
                                           prima.facie.null,
                                           temporal.priority$adj.matrix,temporal.priority$edge.confidence.matrix);

    ## Perform the hypergeometric test for each pair of events.
    
    for (i in 1:ncol(adj.matrix)) {
        for (j in i:nrow(adj.matrix)) {

            ## The diagonal (self cause) and the other invalid edges
            ## have not to be considered.
            
            if (adj.matrix[i, j] != 0
                || adj.matrix[j, i] != 0) {
                
                ## Compute the confidence by hypergeometric test for
                ## both j --> i and i --> j.
                
                probability.raising$edge.confidence.matrix[[3, 1]][i, j] =
                    phyper(joint.probs[i,j] * nrow(dataset),
                           marginal.probs[i] * nrow(dataset),
                           nrow(dataset) - marginal.probs[i] * nrow(dataset),
                           marginal.probs[j] * nrow(dataset),lower.tail = FALSE);

                probability.raising$edge.confidence.matrix[[3, 1]][j, i] =
                    probability.raising$edge.confidence.matrix[[3,1]][i,j];
            } else {
                probability.raising$edge.confidence.matrix[[3, 1]][i, j] = 1;
                probability.raising$edge.confidence.matrix[[3, 1]][j, i] = 1;
            }
        }
    }

    ## Remove any cycle.
    
    adj.matrix.cyclic = probability.raising$adj.matrix
    if (length(temporal.priority$not.ordered) > 0
        || !is.na(hypotheses[1])) {

        if (!silent)
            cat('*** Loop detection found loops to break.\n')

        weights.temporal.priority =
            probability.raising$edge.confidence.matrix[[2, 1]];
        weights.matrix =
            probability.raising$edge.confidence.matrix[[2, 1]] +
                probability.raising$edge.confidence.matrix[[3, 1]];
        acyclic.topology =
            remove.cycles(probability.raising$adj.matrix,
                          weights.temporal.priority,
                          weights.matrix,
                          temporal.priority$not.ordered,
                          hypotheses,
                          silent);
        adj.matrix.acyclic = acyclic.topology$adj.matrix;
    } else {
        adj.matrix.acyclic = probability.raising$adj.matrix;
    }
    adj.matrix =
        list(adj.matrix.cyclic = adj.matrix.cyclic,
             adj.matrix.acyclic=adj.matrix.acyclic)

    ## Save the results and return them.
    
    prima.facie.topology <-
        list(adj.matrix = adj.matrix,
             edge.confidence.matrix = probability.raising$edge.confidence.matrix);
    return(prima.facie.topology);
}


# select the set of the prima facie parents (with bootstrap) for each node
# based on Suppes' definition of causation
# @title get.prima.facie.parents.do.boot
# @param dataset a valid dataset
# @param hypotheses hypotheses object related to dataset
# @param nboot integer number (greater than 0) of bootstrap sampling to be performed
# @param pvalue pvalue for the tests (value between 0 and 1)
# @param adj.matrix adjacency matrix of the initially valid edges
# @param min.boot minimum number of bootstrapping to be performed
# @param min.stat should I keep bootstrapping untill I have nboot valid values?
# @param boot.seed seed to be used for the sampling
# @param silent Should I be verbose?
# @return prima.facie.parents list of the set (if any) of prima facie parents for each node
#
get.prima.facie.parents.do.boot <- function(dataset,
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
        get.bootstapped.scores(dataset,
                               nboot,adj.matrix,
                               min.boot,
                               min.stat,
                               boot.seed,
                               silent);

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
        get.prima.facie.causes.do.boot(adj.matrix,
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


# select the set of the prima facie parents (without bootstrap) for each node based on Suppes' definition of causation
# @title get.prima.facie.parents.no.boot
# @param dataset a valid dataset
# @param hypotheses hypotheses object associated to dataset
# @param adj.matrix adjacency matrix of the initially valid edges
# @param silent Should I be verbose?
# @return prima.facie.parents: list of the set (if any) of prima facie parents for each node
#
get.prima.facie.parents.no.boot <- function(dataset,
                                            hypotheses,
                                            adj.matrix,
                                            silent) {

    ## Compute the scores from the dataset.
    
    scores = get.dag.scores(dataset,adj.matrix);

    ## Remove all the edges not representing a prima facie causes.
    
    prima.facie.topology =
        get.prima.facie.causes.no.boot(adj.matrix,
                                       hypotheses,
                                       scores$marginal.probs,
                                       scores$prima.facie.model,
                                       scores$prima.facie.null,
                                       dataset,
                                       scores$joint.probs,
                                       silent);

    ## Save the results return them.
    
    prima.facie.parents <-
        list(marginal.probs = scores$marginal.probs,
             joint.probs = scores$joint.probs,
             adj.matrix = prima.facie.topology$adj.matrix,
             pf.confidence = prima.facie.topology$edge.confidence.matrix);
    return(prima.facie.parents);
}


# reconstruct the best causal topology by likelihood fit
# @title perform.likelihood.fit.capri
# @param dataset a valid dataset
# @param adj.matrix the adjacency matrix of the prima facie causes
# @param command type of search, either hill climbing (hc) or tabu (tabu)
# @param regularization regularization term to be used in the likelihood fit
# @return topology: the adjacency matrix of both the prima facie and causal topologies
#
perform.likelihood.fit.capri <- function(dataset,
                                         adj.matrix,
                                         command,
                                         regularization) {

    ## Adjacency matrix of the topology reconstructed by likelihood
    ## fit.
    adj.matrix.fit = array(0, c(nrow(adj.matrix), ncol(adj.matrix)))
    rownames(adj.matrix.fit) = colnames(dataset)
    colnames(adj.matrix.fit) = colnames(dataset)

    ## Create a categorical data frame from the dataset.
    data = as.categorical.dataset(dataset)

    # Perform the likelihood fit
    adj.matrix.fit = lregfit(data,
        adj.matrix,
        adj.matrix.fit,
        regularization,
        command)

    ## Save the results and return them.
    
    adj.matrix =
        list(adj.matrix.pf = adj.matrix,
             adj.matrix.fit = adj.matrix.fit);
    topology = list(adj.matrix = adj.matrix);
    return(topology);
}


# remove any cycle from a given cyclic topology
# @title remove.cycles
# @param adj.matrix adjacency matrix of the topology
# @param weights.temporal.priority weighted matrix to be used to remove the cycles involving atomic events
# @param weights.matrix weighted matrix to be used to remove the cycles involving hypotheses
# @param not.ordered list of the nodes to be orderd
# @param hypotheses hypotheses to evaluate potential cycles
# @param silent Should I be verbose?
# @return acyclic.topology: structure representing the best acyclic topology
#
remove.cycles <- function(adj.matrix,
                          weights.temporal.priority,
                          weights.matrix,
                          not.ordered,
                          hypotheses = NA,
                          silent) {

    
    total.edges = length(which(adj.matrix == 1))
    removed = 0

    ## Evaluate the possible cycles involving atomic events.

    if (length(not.ordered) > 0) {

        ## Consider only the edges that were not ordered by temporal
        ## priority.
        
        curr.edge.pos = 0;
        for (i in 1:length(not.ordered)) {

            ## Consider the events i and j.
            
            curr.edge = not.ordered[[i]];
            curr.edge.i = curr.edge[1, 1];
            curr.edge.j = curr.edge[2, 1];

            ## check if i and j still create a cycle.
            
            if (adj.matrix[curr.edge.i, curr.edge.j] == 1
                && adj.matrix[curr.edge.j, curr.edge.i] == 1) {

                ## Get the scores of the two edges.
                
                curr.score.i.j =
                    weights.temporal.priority[curr.edge.i, curr.edge.j];
                curr.score.j.i =
                    weights.temporal.priority[curr.edge.j, curr.edge.i];

                ## Choose an edge based on the score.
                
                if (curr.score.i.j < curr.score.j.i) {
                    
                    ## if i --> j is more confident (lower score) then
                    ## j --> i
                    
                    removed = removed + 1
                    ## cat("Removing edge ",colnames(adj.matrix)[curr.edge.j]," to ",colnames(adj.matrix)[curr.edge.i],"\n");
                    adj.matrix[curr.edge.j, curr.edge.i] = 0;
                } else {
                    ## otherwise
                    removed = removed + 1
                    ## cat("Removing edge ",colnames(adj.matrix)[curr.edge.i]," to ",colnames(adj.matrix)[curr.edge.j],"\n");
                    adj.matrix[curr.edge.i, curr.edge.j] = 0;
                }
            }
        }
    }

    ## Create the structures where to save the weights in increasing
    ## order of confidence.
    
    ordered.weights <- vector();
    ordered.edges <- list();

    ## Consider the patterns related the hypotheses.
    
    if (!is.na(hypotheses[1])) {

        ## If I have hypotheses, add the edges to be evaluated during
        ## the loop removal.
        
        curr.edge.pos = 0;
        for (i in 1:nrow(adj.matrix)) {
            for (j in 1:nrow(adj.matrix)) {
                if (adj.matrix[i, j] == 1) {
                    ordered.weights =
                        rbind(ordered.weights, weights.matrix[i, j]);
                    curr.edge.pos = curr.edge.pos + 1;
                    new.edge <- array(0, c(2, 1));
                    new.edge[1, 1] = i;
                    new.edge[2, 1] = j;
                    ordered.edges[curr.edge.pos] = list(new.edge);
                }
            }
        }

        ## Sort the edges in increasing order of confidence (i.e. the
        ## edges with lower pvalue are the most confident).
        
        ordered.edges =
            ordered.edges[sort(unlist(ordered.weights),
                               decreasing = TRUE,
                               index.return = TRUE)$ix];
    }

    ## Visit the ordered edges and remove the ones that are causing
    ## any cycle.
    
    if (length(ordered.edges) > 0) {

        ## Expanded matrix to be considered in removing the loops.
        
        expansion =
            hypotheses.expansion(input_matrix = adj.matrix,
                                 map=hypotheses$hstructure,
                                 expand = TRUE,
                                 skip.disconnected = FALSE);


        for (i in 1:length(ordered.edges)) {

            ## Consider the edge i-->j
            curr.edge = ordered.edges[[i]];
            curr.edge.i = curr.edge[1,1];
            curr.edge.j = curr.edge[2,1];

            ## Resolve the mapping from the adj.matrix to the expanded
            ## one both for curr.edge.i and curr.edge.j
            
            if (colnames(adj.matrix)[curr.edge.i] %in% expansion[[2]]) {
                curr.edge.i.exp =
                    which(colnames(expansion[[1]]) %in%
                          names(expansion[[2]])[which(expansion[[2]] %in%
                                                      colnames(adj.matrix)[curr.edge.i])]);
            } else {
                curr.edge.i.exp =
                    which(colnames(expansion[[1]]) %in%
                          colnames(adj.matrix)[curr.edge.i]);
            }
            
            if (colnames(adj.matrix)[curr.edge.j] %in% expansion[[2]]) {
                curr.edge.j.exp =
                    which(colnames(expansion[[1]]) %in%
                          names(expansion[[2]])[which(expansion[[2]] %in%
                                                      colnames(adj.matrix)[curr.edge.j])]);
            } else {
                curr.edge.j.exp =
                    which(colnames(expansion[[1]]) %in%
                          colnames(adj.matrix)[curr.edge.j]);
            }

            ## Search for loops between curr.edge.i and curr.edge.j
            
            curr.graph = graph.adjacency(expansion[[1]], mode = "directed")

            is.path = suppressWarnings(get.shortest.paths(curr.graph,
                                       curr.edge.j.exp,
                                       curr.edge.i.exp)$vpath)

            is.path = length(unlist(is.path))

            ## If there is a path between the two nodes, remove edge i --> j
            
            if (is.path > 0) {
                removed = removed + 1

                ## cat("Removing edge ",colnames(adj.matrix)[curr.edge.i]," to ",colnames(adj.matrix)[curr.edge.j],"\n");

                expansion[[1]][curr.edge.i.exp, curr.edge.j.exp] = 0;
                adj.matrix[curr.edge.i, curr.edge.j] = 0;
            }
        }

        if (!silent)
            cat(paste0('\tRemoved ',
                       removed,
                       ' edges out of ',
                       total.edges,
                       ' (',
                       round(100 * removed/total.edges, 0),
                       '%)\n'))
    }

    ## Save the results and return them.
    
    acyclic.topology = list(adj.matrix = adj.matrix);
    return(acyclic.topology);
}


# verify the probability raising condition
# @title verify.probability.raising.do.boot
# @param prima.facie.model.distributions distributions of the prima facie model
# @param prima.facie.null.distributions distributions of the prima facie null
# @param pvalue minimum pvalue for the Mann-Whitney U tests to be significant
# @param adj.matrix adjacency matrix of the topology
# @param edge.confidence.matrix matrix of the confidence of each edge
# @return probability.raising: list describing the causes where probability raising is verified
#
verify.probability.raising.do.boot <- function(prima.facie.model.distributions,
                                               prima.facie.null.distributions,
                                               pvalue,
                                               adj.matrix,
                                               edge.confidence.matrix) {

    ## Evaluate the probability raising condition.
    
    for (i in 1:nrow(adj.matrix)) {
        for (j in i:ncol(adj.matrix)) {

            ## The diagonal (self cause) and the other invalid edges
            ## have not to be considered.
            
            if (adj.matrix[i, j] != 0
                || adj.matrix[j, i] != 0) {

                ## pvalue for the probability raising condition for i --> j

                second.pvalue.i.j = suppressWarnings(
                    wilcox.test(unlist(prima.facie.model.distributions[i, j]),
                                unlist(prima.facie.null.distributions[i,j]),
                                alternative = "greater",
                                mu = 0))$p.value;
                if (is.na(second.pvalue.i.j)
                    || is.nan(second.pvalue.i.j)) {

                    ## In this case the two distributions are exactly
                    ## identical.
                    
                    second.pvalue.i.j = 1;
                }

                ## In this case i --> j is not valid
                
                if (second.pvalue.i.j >= pvalue) {
                    adj.matrix[i, j] = 0;
                }

                ## pvalue for the probability raising condition for j --> i

                second.pvalue.j.i = suppressWarnings(
                    wilcox.test(unlist(prima.facie.model.distributions[j, i]),
                                unlist(prima.facie.null.distributions[j,i]),
                                alternative = "greater", mu = 0))$p.value;
                if (is.na(second.pvalue.j.i)
                    || is.nan(second.pvalue.j.i)) {

                    ## In this case the two distributions are exactly
                    ## identical.
                    
                    second.pvalue.j.i = 1;
                }

                ## In this case j --> i is not valid
                
                if (second.pvalue.j.i >= pvalue) {
                    adj.matrix[j, i] = 0;
                }

                ## Save the confidence for i--> j and j --> i

                tmp = edge.confidence.matrix[[2, 1]];
                tmp[i, j] = second.pvalue.i.j;
                tmp[j, i] = second.pvalue.j.i;
                edge.confidence.matrix[2, 1] = list(tmp);
            } else {
                tmp = edge.confidence.matrix[[2, 1]];
                tmp[i, j] = 1;
                tmp[j, i] = 1;
                edge.confidence.matrix[2, 1] = list(tmp);
            }
        }
    }

    ## Save the results and return them.
    
    probability.raising <-
        list(adj.matrix = adj.matrix,
             edge.confidence.matrix = edge.confidence.matrix);
    return(probability.raising);
}


# verify the probability raising condition without bootstrap
# @title verify.probability.raising.no.boot
# @param prima.facie.model prima facie model
# @param prima.facie.null prima facie null
# @param adj.matrix adjacency matrix of the topology
# @param edge.confidence.matrix matrix of the confidence of each edge
# @return probability.raising: adjacency matrix where temporal priority is verified

verify.probability.raising.no.boot <- function(prima.facie.model,
                                               prima.facie.null,
                                               adj.matrix,
                                               edge.confidence.matrix) {

    ## Evaluate the probability raising condition.
    
    for (i in 1:nrow(adj.matrix)) {
        for (j in i:ncol(adj.matrix)) {

            ## The diagonal (self cause) and the other invalid edges
            ## have not to be considered probability raising
            ## condition: if P(j|i)>P(j|not i) the edge i --> j is
            ## valid for probability raising.
            
            if (adj.matrix[i, j] != 0
                || adj.matrix[j, i] != 0) {

                ## In this case i --> j is not valid
                
                if (prima.facie.model[i, j] <= prima.facie.null[i, j]) {
                    adj.matrix[i, j] = 0;
                }

                ## In this case j --> i is not valid
                if (prima.facie.model[j, i] <= prima.facie.null[j, i]) {
                    adj.matrix[j, i] = 0;
                }

                ## Save the confidence for i-->j and j --> i

                tmp = edge.confidence.matrix[[2, 1]];
                tmp[i, j] = min(prima.facie.null[i, j] / prima.facie.model[i, j],1);
                tmp[j, i] = min(prima.facie.null[j, i] / prima.facie.model[j, i],1);
                edge.confidence.matrix[2, 1] = list(tmp);
            } else {
                tmp = edge.confidence.matrix[[2, 1]];
                tmp[i, j] = 1;
                tmp[j, i] = 1;
                edge.confidence.matrix[2, 1] = list(tmp);
            }
        }
    }

    ## Save the results and return them.
    
    probability.raising <-
        list(adj.matrix = adj.matrix,
             edge.confidence.matrix = edge.confidence.matrix);
    return(probability.raising);
}


# verify the temporal priority condition with bootstrap
# @title verify.temporal.priority.do.boot
# @param marginal.probs.distributions distributions of the bootstrapped marginal probabilities
# @param pvalue minimum pvalue for the Mann-Whitney U tests to be significant
# @param adj.matrix adjacency matrix of the topology
# @param edge.confidence.matrix matrix of the confidence of each edge
# @return temporal.priority: list describing the causes where temporal priority is verified
#
verify.temporal.priority.do.boot <- function(marginal.probs.distributions,
                                             pvalue,
                                             adj.matrix,
                                             edge.confidence.matrix) {

    ## Evalutate the temporal priority condition for each pair of
    ## edges.
    
    not.ordered = list();
    counter = 0;
    for (i in 1:nrow(adj.matrix)) {
        for (j in i:ncol(adj.matrix)) {

            ## The diagonal (self cause) and the other invalid edges
            ## have not to be considered.
            
            if (adj.matrix[i, j] != 0
                || adj.matrix[j, i] != 0) {

                ## [i,j] refers to causation i --> j
                ## Temporal priority condition: if P(i) > P(j) the edge
                ## i --> j is valid for temporal priority.
                
                ## Test i --> j
                first.pvalue.i.j = suppressWarnings(
                    wilcox.test(unlist(marginal.probs.distributions[i, 1]),
                                unlist(marginal.probs.distributions[j, 1]),
                                alternative = "greater",
                                mu = 0))$p.value;
                if (is.na(first.pvalue.i.j)
                    || is.nan(first.pvalue.i.j)) {

                    ## In this case the two distributions are exactly
                    ## identical.
                    
                    first.pvalue.i.j = 1;
                }

                ## Test j --> i
                
                first.pvalue.j.i = suppressWarnings(
                    wilcox.test(unlist(marginal.probs.distributions[j, 1]),
                                unlist(marginal.probs.distributions[i, 1]),
                                alternative = "greater",
                                mu = 0))$p.value;
                if (is.na(first.pvalue.j.i)
                    || is.nan(first.pvalue.j.i)) {
                    
                    ## In this case the two distributions are exactly
                    ## identical.
                    
                    first.pvalue.j.i = 1;
                }

                ## In this case i is before j and j --> i is not valid.

                if (first.pvalue.j.i >= pvalue
                    && first.pvalue.i.j < pvalue) {

                    ## [j,i] = 0 means j is after i, i.e. it can not be causing i
                    adj.matrix[j,i] = 0;
                }

                ## In this case j is before i and i --> j is not valid.

                else if (first.pvalue.j.i < pvalue
                         && first.pvalue.i.j >= pvalue) {
                    ## [i,j] = 0 means i is after j, i.e. it can not be causing j
                    adj.matrix[i,j] = 0;
                }
                
                ## In this case, a total time order between i and j
                ## can not be defined.
                
                else {
                    
                    ## No temporal priority induced by the topology
                    ## can be inferred.
                    
                    counter = counter + 1;
                    curr.not.ordered = array(-1, c(2, 1));
                    curr.not.ordered[1, 1] = i;
                    curr.not.ordered[2, 1] = j;
                    not.ordered[counter] = list(curr.not.ordered);
                }

                ## Save the confidence for i --> j and j --> i

                tmp = edge.confidence.matrix[[1, 1]];
                tmp[i, j] = first.pvalue.i.j;
                tmp[j, i] = first.pvalue.j.i;
                edge.confidence.matrix[1, 1] = list(tmp);
            } else {
                tmp = edge.confidence.matrix[[1, 1]];
                tmp[i, j] = 1;
                tmp[j, i] = 1;
                edge.confidence.matrix[1, 1] = list(tmp);
            }
        }
    }

    ## Save the results and return them.
    
    temporal.priority <-
        list(adj.matrix = adj.matrix,
             edge.confidence.matrix = edge.confidence.matrix,
             not.ordered = not.ordered);
    return(temporal.priority);
}


# verify the temporal priority condition without bootstrap
# @title verify.temporal.priority.no.boot
# @param marginal.probs marginal probabilities
# @param adj.matrix adjacency matrix of the topology
# @param edge.confidence.matrix matrix of the confidence of each edge
# @return temporal.priority: adjacency matrix where temporal priority is verified
#
verify.temporal.priority.no.boot <- function(marginal.probs,
                                             adj.matrix,
                                             edge.confidence.matrix) {

    ## Evalutate the temporal priority condition for each pair of
    ## edges.
    
    not.ordered = list();
    counter = 0;
    for (i in 1:nrow(adj.matrix)) {
        for (j in i:ncol(adj.matrix)) {

            ## The diagonal (self cause) and the other invalid edges
            ## have not to be considered.
            
            if (adj.matrix[i,j] != 0
                || adj.matrix[j, i] != 0) {

                ## [i,j] refers to causation i --> j

                ## Temporal priority condition: if P(i)>P(j) the edge
                ## i --> j is valid for temporal priority in this case
                ## i is before j and j --> i is not valid.
                
                if (marginal.probs[i, 1] > marginal.probs[j, 1]) {
                    ## [j,i] = 0 means j is after i, i.e. it can not be causing i

                    adj.matrix[j, i] = 0;
                }

                ## In this case j is before i and i --> j is not valid
                
                else if (marginal.probs[j, 1] > marginal.probs[i, 1]) {
                    ## [i,j] = 0 means i is after j, i.e. it can not be causing j

                    adj.matrix[i,j] = 0;
                }
                
                ## In this case, a total time order between i and j
                ## can not be defined.
                
                else {
                    
                    ## No temporal priority induced by the topology
                    ## can be inferred.
                    
                    counter = counter + 1;
                    curr.not.ordered = array(-1, c(2, 1));
                    curr.not.ordered[1, 1] = i;
                    curr.not.ordered[2, 1] = j;
                    not.ordered[counter] = list(curr.not.ordered);
                }

                ## Save the confidence for i --> j and j --> i

                tmp = edge.confidence.matrix[[1, 1]];
                tmp[i, j] = min(marginal.probs[j, 1] / marginal.probs[i, 1], 1);
                tmp[j, i] = min(marginal.probs[i,1]/marginal.probs[j,1],1);
                edge.confidence.matrix[1, 1] = list(tmp);

            } else {
                tmp = edge.confidence.matrix[[1, 1]];
                tmp[i, j] = 1;
                tmp[j, i] = 1;
                edge.confidence.matrix[1, 1] = list(tmp);
            }
        }
    }

    ## Save the results and return them.
    
    temporal.priority <-
        list(adj.matrix = adj.matrix,
             edge.confidence.matrix = edge.confidence.matrix,
             not.ordered = not.ordered);
    return(temporal.priority);
}


#### end of file -- capri.algorithm.R
