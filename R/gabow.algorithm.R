#### TRONCO: a tool for TRanslational ONCOlogy
####
#### Copyright (c) 2015-2016, Marco Antoniotti, Giulio Caravagna, Luca De Sano,
#### Alex Graudenzi, Giancarlo Mauri, Bud Mishra and Daniele Ramazzotti.
####
#### All rights reserved. This program and the accompanying materials
#### are made available under the terms of the GNU GPL v3.0
#### which accompanies this distribution.

# reconstruct the best topology based on probabilistic causation and maximum likelihood estimation
# @title gabow.fit
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
gabow.fit <- function(dataset,
                        regularization = "no_reg",
                        score = "pmi",
                        do.boot = TRUE,
                        nboot = 100,
                        pvalue = 0.05,
                        min.boot = 3,
                        min.stat = TRUE,
                        boot.seed = NULL,
                        silent = FALSE,
                        epos = 0.0,
                        eneg = 0.0,
                        do.raising = FALSE ) {

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
    
    valid.dataset = check.dataset(dataset, adj.matrix, FALSE, epos, eneg)
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
        
    # set the input prior matrix for Gabow
    if(do.raising==TRUE) {
        adj.matrix.gabow = prima.facie.parents$adj.matrix$adj.matrix.cyclic
    }
    else {
        adj.matrix.gabow = prima.facie.parents$adj.matrix$adj.matrix.cyclic.tp
    }

    ## Perform the likelihood fit with the required strategy.
    
    model = list();
    for (reg in regularization) {
        for(my_score in score) {

            ## Perform the likelihood fit with the chosen regularization
            ## score on the prima facie topology.
            
            if (!silent)
                cat('*** Performing likelihood-fit with regularization:', reg, '.\n')
            best.parents =
                perform.likelihood.fit.gabow(dataset,
                                       adj.matrix.gabow,
                                       regularization = reg,
                                       score = my_score,
                                       marginal.probs = prima.facie.parents$marginal.probs,
                                       joint.probs = prima.facie.parents$joint.probs)
    
            ## Set the structure to save the conditional probabilities of
            ## the reconstructed topology.
    
            reconstructed.model = create.model(dataset,
                best.parents,
                prima.facie.parents)
    
            model.name = paste('gabow', reg, my_score, sep='_')
            model[[model.name]] = reconstructed.model
            
        }
    }

    ## Set the execution parameters.
    
    parameters =
        list(algorithm = "GABOW",
             regularization = regularization,
             score = score,
             do.boot = do.boot,
             nboot = nboot,
             pvalue = pvalue,
             min.boot = min.boot,
             min.stat = min.stat,
             boot.seed = boot.seed,
             silent = silent,
             error.rates = list(epos=epos,eneg=eneg),
             do.raising = do.raising)

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
# @title perform.likelihood.fit.gabow
# @param dataset a valid dataset
# @param adj.matrix the adjacency matrix of the prima facie causes
# @param regularization regularization term to be used in the likelihood fit
# @param command type of search, either hill climbing (hc) or tabu (tabu)
# @return topology: the adjacency matrix of both the prima facie and causal topologies
#
perform.likelihood.fit.gabow = function(dataset,
                                        adj.matrix,
                                        regularization,
                                        command = "hc",
                                        marginal.probs,
                                        joint.probs) {

    data = as.categorical.dataset(dataset)
    adj.matrix.prima.facie = adj.matrix
    
    # adjacency matrix of the topology reconstructed by likelihood fit
    adj.matrix.fit = array(0,c(nrow(adj.matrix),ncol(adj.matrix)))
    rownames(adj.matrix.fit) = colnames(dataset)
    colnames(adj.matrix.fit) = colnames(dataset)
    
    # set at most one parent per node based on the given score and Gabow search
    
    # 
    curr_nodes = which(apply(adj.matrix,1,sum)==0)
    curr_adj.matrix.fit = adj.matrix.fit
    visited_nodes = as.vector(curr_nodes)
    while(length(curr_nodes)>0) {
        
        new_curr_nodes = NULL
        
        # visit the curr_nodes
        for (i in curr_nodes) {
            
            # get the candidate parents of i
            curr_parents = as.vector(which(adj.matrix[,i] == 1))
            
            # if I have more then one valid parent
            if (length(curr_parents) > 0) {
                
                # find the best parent
                curr_best_parent = -1
                curr_best_score = NA
                
                for (j in curr_parents) {
                    
                    # consider the network augmenting the current one with the edge j --> i
                    tmp_adj.matrix.fit =  curr_adj.matrix.fit
                    tmp_adj.matrix.fit[j,i] = 1
                    colnames(tmp_adj.matrix.fit) = colnames(data)
                    rownames(tmp_adj.matrix.fit) = colnames(data)
                    
                    # create a bnlearn object from tmp_adj.matrix.fit
                    curr.bayes.net = NULL
                    curr.bayes.net$data = data
                    net = empty.graph(colnames(tmp_adj.matrix.fit))
                    for (from in rownames(tmp_adj.matrix.fit)) {
                        for (to in colnames(tmp_adj.matrix.fit)) {
                            if (tmp_adj.matrix.fit[from, to] == 1) {
                                net = set.arc(net, from, to)
                            }
                        }
                    }
                    curr.bayes.net$net = net
                    
                    # evaluate any candidate parent in terms of likelihood
                    new_score = logLik(curr.bayes.net$net, data = curr.bayes.net$data)
                    
                    # compare the likelihood of the current candidate parent 
                    # with the current best one
                    if (is.na(curr_best_score) || (new_score > curr_best_score)) {
                        curr_best_parent = j
                        curr_best_score = new_score
                    }
                
                }
                
                # set the best parent in the inferred matrix
                curr_adj.matrix.fit[curr_best_parent,i] = 1
                
                # set the current best parent in the new list of curr nodes
                new_curr_nodes = c(new_curr_nodes,curr_best_parent)
                
            }
            
        }
        
        # find out the nodes to still be visited within new_curr_nodes
        new_curr_nodes = unique(new_curr_nodes)
        new_idx = which(!(new_curr_nodes%in%visited_nodes))
        if(length(new_idx)>0) {
            new_curr_nodes = new_curr_nodes[new_idx]
        }
        else {
            new_curr_nodes = NULL
        }
        
        # set the list of the next curr_nodes
        curr_nodes = new_curr_nodes
        
        # update the list of visited nodes
        visited_nodes = c(visited_nodes,curr_nodes)
        
    }
    
    # consider any not visited node
    not_visited = which(!((1:nrow(adj.matrix))%in%visited_nodes))
    
    if(length(not_visited)>0) {
        
        for (i in not_visited) {
            
            # get the candidate parents of i
            curr_parents = as.vector(which(adj.matrix[,i] == 1))
            
            # if I have more then one valid parent
            if (length(curr_parents) > 0) {
                
                # find the best parent
                curr_best_parent = -1
                curr_best_score = NA
                
                for (j in curr_parents) {
                    
                    # consider the network augmenting the current one with the edge j --> i
                    tmp_adj.matrix.fit =  curr_adj.matrix.fit
                    tmp_adj.matrix.fit[j,i] = 1
                    colnames(tmp_adj.matrix.fit) = colnames(data)
                    rownames(tmp_adj.matrix.fit) = colnames(data)
                    
                    # create a bnlearn object from tmp_adj.matrix.fit
                    curr.bayes.net = NULL
                    curr.bayes.net$data = data
                    net = empty.graph(colnames(tmp_adj.matrix.fit))
                    for (from in rownames(tmp_adj.matrix.fit)) {
                        for (to in colnames(tmp_adj.matrix.fit)) {
                            if (tmp_adj.matrix.fit[from, to] == 1) {
                                net = set.arc(net, from, to)
                            }
                        }
                    }
                    curr.bayes.net$net = net
                    
                    # evaluate any candidate parent in terms of likelihood
                    new_score = logLik(curr.bayes.net$net, data = curr.bayes.net$data)
                    
                    # compare the likelihood of the current candidate parent 
                    # with the current best one
                    if (is.na(curr_best_score) || (new_score > curr_best_score)) {
                        curr_best_parent = j
                        curr_best_score = new_score
                    }
                
                }
                
                # set the best parent in the inferred matrix
                curr_adj.matrix.fit[curr_best_parent,i] = 1
                
            }
            
        }
        
    }
    
    # perform the likelihood fit if requested
    adj.matrix.fit = lregfit(data,
        curr_adj.matrix.fit,
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


#### end of file -- gabow.algorithm.R
