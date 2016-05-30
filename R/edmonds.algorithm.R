#### TRONCO: a tool for TRanslational ONCOlogy
####
#### Copyright (c) 2015-2016, Marco Antoniotti, Giulio Caravagna, Luca De Sano,
#### Alex Graudenzi, Giancarlo Mauri, Bud Mishra and Daniele Ramazzotti.
####
#### All rights reserved. This program and the accompanying materials
#### are made available under the terms of the GNU GPL v3.0
#### which accompanies this distribution.

# reconstruct the best topology based on probabilistic causation and Edmonds algorithm
# @title edmonds.fit
# @param dataset a dataset describing a progressive phenomenon
# @param regularization regularizators to be used for the likelihood fit
# @param score the score to be used, could be either pointwise mutual information (pmi) or conditional entropy (entropy)
# @param do.boot should I perform bootstrap? Yes if TRUE, no otherwise
# @param nboot integer number (greater than 0) of bootstrap sampling to be performed
# @param pvalue pvalue for the tests (value between 0 and 1)
# @param min.boot minimum number of bootstrapping to be performed
# @param min.stat should I keep bootstrapping untill I have nboot valid values?
# @param boot.seed seed to be used for the sampling
# @param silent should I be verbose?
# @param epos error rate of false positive errors
# @param eneg error rate of false negative errors
# @return topology: the reconstructed tree topology
#
edmonds.fit <- function(dataset,
                        regularization = "no_reg",
                        score = "pmi",
                        do.boot = TRUE,
                        nboot = 100,
                        pvalue = 0.05,
                        min.boot = 3,
                        min.stat = TRUE,
                        boot.seed = NULL,
                        silent = FALSE,
                        epos,
                        eneg ) {

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
                                            silent,
                                            epos,
                                            eneg);
    } else {
        if (!silent)
            cat('*** Computing selective advantage scores (prima facie).\n')
        prima.facie.parents =
            get.prima.facie.parents.no.boot(dataset,
                                            NA,
                                            adj.matrix,
                                            silent,
                                            epos,
                                            eneg);
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
        for(my_score in score) {

            ## Perform the likelihood fit with the chosen regularization
            ## and the chosed score on the prima facie topology.
            
            if (!silent)
                cat('*** Performing likelihood-fit with regularization:', reg, 'and score:', my_score, '.\n')
            best.parents =
                perform.likelihood.fit.edmonds(dataset,
                                       adj.matrix.prima.facie,
                                       regularization = reg,
                                       score = my_score,
                                       marginal.probs = prima.facie.parents$marginal.probs,
                                       joint.probs = prima.facie.parents$joint.probs)
    
            ## Set the structure to save the conditional probabilities of
            ## the reconstructed topology.
    
            reconstructed.model = create.model(dataset,
                best.parents,
                prima.facie.parents)
    
            model.name = paste('edmonds', reg, my_score, sep='_')
            model[[model.name]] = reconstructed.model
            
        }
    }

    ## Set the execution parameters.
    
    parameters =
        list(algorithm = "EDMONDS",
             regularization = regularization,
             score = score,
             do.boot = do.boot,
             nboot = nboot,
             pvalue = pvalue,
             min.boot = min.boot,
             min.stat = min.stat,
             boot.seed = boot.seed,
             silent = silent,
             error.rates = list(epos=epos,eneg=eneg));

    ## Return the results.
    
    topology =
        list(dataset = dataset,
             hypotheses = NA,
             adj.matrix.prima.facie = adj.matrix.prima.facie,
             adj.matrix.prima.facie.cyclic = prima.facie.parents$adj.matrix$adj.matrix.cyclic,
             confidence = prima.facie.parents$pf.confidence,
             model = model,
             parameters = parameters,
             execution.time = (proc.time() - ptm))
    topology = rename.reconstruction.fields(topology, dataset)
    return(topology)
}


# reconstruct the best causal topology by Edmonds algorithm combined with probabilistic causation
# @title perform.likelihood.fit.edmonds
# @param dataset a valid dataset
# @param adj.matrix the adjacency matrix of the prima facie causes
# @param regularization regularization term to be used in the likelihood fit
# @param score the score to be used by edmonds algorithm. Could be either pointwise mutual information (pmi) or conditional entropy (entropy)
# @param command type of search, either hill climbing (hc) or tabu (tabu)
# @return topology: the adjacency matrix of both the prima facie and causal topologies
#
perform.likelihood.fit.edmonds = function(dataset,
                                          adj.matrix,
                                          regularization,
                                          score,
                                          command = "hc",
                                          marginal.probs,
                                          joint.probs){

    data = as.categorical.dataset(dataset)
    adj.matrix.prima.facie = adj.matrix
    
    # adjacency matrix of the topology reconstructed by likelihood fit
    adj.matrix.fit = array(0,c(nrow(adj.matrix),ncol(adj.matrix)))
    rownames(adj.matrix.fit) = colnames(dataset)
    colnames(adj.matrix.fit) = colnames(dataset)
       
    # set at most one parent per node based on mutual information
    for (i in 1:ncol(adj.matrix)) {
            
        # consider the parents of i
        curr_parents = which(adj.matrix[,i] == 1)
        
        # if I have more then one valid parent
        if (length(curr_parents) > 1) {
            
            # find the best parent
            curr_best_parent = -1
            curr_best_score = -1
            for (j in curr_parents) {
                
                # if the event is valid
                if(joint.probs[i,j]>=0) {
                    # compute the chosen score
                    new_score = compute.edmonds.score(joint.probs[i,j],marginal.probs[i],marginal.probs[j],score)
                }
                # else, if the two events are indistinguishable
                else if(joint.probs[i,j]<0) {
                    new_score = Inf
                }
                
                if (new_score > curr_best_score) {
                    curr_best_parent = j
                    curr_best_score = new_score
                }
            }
            
            # set the best parent
            for (j in curr_parents) {
                if (j != curr_best_parent) {
                    adj.matrix[j,i] = 0
                }
            }
        }
    }
    
    # perform the likelihood fit if requested
    adj.matrix.fit = lregfit(data,
        adj.matrix,
        adj.matrix.fit,
        regularization,
        command)
    
    ## Save the results and return them.
    
    adj.matrix =
        list(adj.matrix.pf = adj.matrix.prima.facie,
             adj.matrix.fit = adj.matrix.fit)
    topology = list(adj.matrix = adj.matrix)
    return(topology)

}

# compute either pointwise mutual information (pmi) or conditional entropy (mle) for the edge j --> i
compute.edmonds.score = function( joint.prob.i.j, marginal.prob.i, marginal.prob.j, score ) {
    
    # variable to save the score
    new_score = NULL
    
    # this is the pointwise mutual information
    if(score=="pmi") {
        
        # compute the pointwise mutual information for i and j
        # that is log(P(i,j)/[P(i)*P(j)])
        new_score = log(joint.prob.i.j/(marginal.prob.i*marginal.prob.j))
        if(is.nan(new_score)) {
            new_score = 0
        }
        
    }
    # this is the conditional entropy of i given j
    else if(score=="entropy") {
        
        # compute the 4 components of the conditional entropy
        h.i.j = joint.prob.i.j * 
                    log(marginal.prob.j/joint.prob.i.j)
        if(is.nan(h.i.j)) {
            h.i.j = 0
        }
        h.i.not.j = (marginal.prob.i - joint.prob.i.j)  * 
                    log((1-marginal.prob.j)/(marginal.prob.i - joint.prob.i.j))
        if(is.nan(h.i.not.j)) {
            h.i.not.j = 0
        }
        h.not.i.j = (marginal.prob.j - joint.prob.i.j) * 
                    log(marginal.prob.j/(marginal.prob.j - joint.prob.i.j))
        if(is.nan(h.not.i.j)) {
            h.not.i.j = 0
        }
        h.not.i.not.j = (1 - marginal.prob.i - marginal.prob.j + joint.prob.i.j) * 
                    log((1-marginal.prob.j)/(1 - marginal.prob.i - marginal.prob.j + joint.prob.i.j))
        if(is.nan(h.not.i.not.j)) {
            h.not.i.not.j = 0
        }
        
        # compute the entropy
        new_score = h.i.j + h.i.not.j + h.not.i.j + h.not.i.not.j
        
        # this is a maximization problem, hence we use the negate of the conditional entropy
        new_score = - new_score
        
    }
    # this is the pointwise mutual information of i and j adjusted for the normalized 
    # pointwise mutual information of i and not j
    else if(score=="cpmi") {
        
        # compute the first part of the score,
        # i.e., the pointwise mutual information for i and j
        # that is log(P(i,j)/[P(i)*P(j)])
        new_score = log(joint.prob.i.j/(marginal.prob.i*marginal.prob.j))
        if(is.nan(new_score)) {
            new_score = 0
        }
        
        # now compute the second part,
        # i.e., the normalized pointwise mutual information for i and not j
        if(new_score>0) {
            
            # compute first the pointwise mutual information for i and not j
            # that is log(P(i,not j)/[P(i)*P(not j)])
            pmi_i_not_j = log((marginal.prob.i-joint.prob.i.j)/(marginal.prob.i*(1-marginal.prob.j)))
            
            # compute the normalization factor for pmi_i_not_j,
            # that is -log(P(i,not j))
            norm_pmi_i_not_j = - log(marginal.prob.i-joint.prob.i.j)
            
            # compute the normalized pointwise mutual information for i and not j
            npmi_i_not_j = pmi_i_not_j/norm_pmi_i_not_j
            # now I correct for any NA (e.g., -Inf/Inf)
            if(is.nan(npmi_i_not_j)) {
                npmi_i_not_j = -1
            }
            
            # now I can correct new_score for npmi_i_not_j
            new_score = new_score * (- npmi_i_not_j)
            
        }
        
    }
    
    return(new_score)
    
}


#### end of file -- edmonds.algorithm.R
